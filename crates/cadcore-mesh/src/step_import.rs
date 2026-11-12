//! Reading a STEP file that some other CAD system wrote.
//!
//! The older [`crate::step_reader`] + [`crate::tessellate`] pair was built to
//! re-open our *own* exports. It works by **recognising shapes**: it looks at a
//! face, decides "this is a disk", "this is a fillet torus", "this is a full
//! cylinder", and then draws that idealised primitive. Every recognition rule
//! encodes an assumption our writer happens to satisfy — a fillet whose bounding
//! circles have exactly the minor radius, a cylinder capped by two circles, a
//! planar face whose outline is convex.
//!
//! A STEP from anywhere else breaks those assumptions immediately, and the
//! failure mode is the worst kind: it does not error, it *draws something*. A
//! toroidal face whose two rims are not recognised becomes a **complete donut**.
//! A cylinder whose rims are ellipses rather than circles becomes a **10 mm long
//! full tube** (that number is literally hard-coded in the fallback). Arcs
//! become whole circles, because the reader never reads the trim parameters at
//! all and stamps `t_end = 2π` on every circular edge. The result looks like a
//! corrupted model, and the file is fine.
//!
//! This module does the thing a CAD kernel actually does: **trimmed-surface
//! tessellation**. No face is ever recognised as anything. For each face we
//!
//! 1. take the carrier surface and its own `(u, v)` parametrisation,
//! 2. walk the real edge loops, evaluating each curve **between its real trim
//!    parameters** — the arc that ends at the two vertices, not the full circle,
//! 3. map that boundary into `(u, v)`, unwrapping the seam for periodic
//!    surfaces,
//! 4. triangulate the region the boundary actually encloses,
//! 5. map back to 3-D through the surface, so curvature is exact rather than
//!    inferred.
//!
//! What comes out is the face the file describes, whatever shape that is.
//!
//! ## Watertightness
//!
//! Two faces that share an edge must be given **the same points** along it, or
//! the model shows cracks. So every `EDGE_CURVE` is discretised exactly once and
//! cached; the second face to use it gets the same polyline reversed. Interior
//! refinement never moves a boundary vertex, and where it subdivides a boundary
//! edge it splits the *chord*, not the curve, so the shared polyline is
//! unchanged. The importer then counts the mesh's own open edges and reports the
//! number, so "watertight" is a measurement rather than a hope.
//!
//! ## What is not guessed
//!
//! A face whose boundary the importer cannot make sense of is **not drawn** and
//! is counted in [`ImportReport::faces_failed`] with the reason. Silently
//! emitting a plausible-looking primitive instead is exactly the bug this module
//! replaces.

use std::collections::{HashMap, HashSet};

use cadcore_math::{Frame3, Point3, UnitVec3, Vec3};
use serde::Serialize;

use crate::tessellate::FaceMesh;

thread_local! {
    /// Entity id of the face being meshed, so a dump can name it.
    static DUMP_FACE: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
    /// Boundary loops discarded because they lie outside the face they were
    /// attached to. Counted rather than hidden: it is a fault in the file.
    static STRAY_LOOPS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

const TAU: f64 = std::f64::consts::TAU;
const PI: f64 = std::f64::consts::PI;

// ─────────────────────────────── public API ────────────────────────────────

/// Knobs for the import. The defaults suit viewing a scaffold at 0.05 mm.
#[derive(Clone, Debug)]
pub struct ImportOptions {
    /// Maximum distance, in millimetres, between the mesh and the true surface.
    /// Drives both edge discretisation and interior refinement.
    pub chord_tol: f64,
    /// Ceiling on triangles produced for a single face. Refinement stops here
    /// even if the tolerance is not met; the shortfall is reported.
    pub max_triangles_per_face: usize,
}

impl Default for ImportOptions {
    fn default() -> Self {
        Self {
            chord_tol: 0.02,
            max_triangles_per_face: 20_000,
        }
    }
}

/// What the import actually managed to do. Every field is a count taken from
/// the result, not an intention.
#[derive(Clone, Debug, Default, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ImportReport {
    pub solids: usize,
    pub faces_total: usize,
    pub faces_meshed: usize,
    pub faces_failed: usize,
    pub triangles: usize,
    pub vertices: usize,
    /// Millimetres per file unit, read from the file's own unit context.
    pub length_scale_mm: f64,
    /// Triangle edges used by exactly one triangle. Zero for a closed solid;
    /// anything else means the mesh has holes, and how many.
    pub open_edges: usize,
    /// Edges whose two triangles run along them the same way instead of
    /// opposite ways. Zero for a consistently oriented solid; anything else
    /// means part of the model is inside-out.
    pub inconsistent_edges: usize,
    /// Openings left by loops the file put on the wrong face, closed with a
    /// patch so the mesh is a solid again. Each is smaller than a filament and
    /// lies inside the material, so a section never draws one.
    pub capped_holes: usize,
    /// Triangles whose winding disagreed with the surface normal the file
    /// states, and which were turned round to match. A face triangulated in its
    /// own parameter space comes out wound whichever way that space happens to
    /// run; the normal is the file's own statement of which side is out, so it
    /// is the one to believe.
    pub rewound_triangles: usize,
    /// Boundary loops the file attached to a face they are not on, and which
    /// were therefore ignored. Non-zero means the STEP is malformed, not that
    /// the import went wrong.
    pub stray_loops: usize,
    /// Faces the importer refused to draw, grouped by reason with counts.
    pub failures: Vec<String>,
    /// Anything else worth knowing — unsupported entities, ignored assemblies.
    pub warnings: Vec<String>,
    pub bbox_min: [f64; 3],
    pub bbox_max: [f64; 3],
}

/// Read a STEP file and tessellate every face it contains.
///
/// Lengths come out in millimetres regardless of the file's own unit.
pub fn import_step(contents: &str, opts: &ImportOptions) -> Result<(Vec<FaceMesh>, ImportReport), String> {
    let doc = Doc::parse(contents);
    if doc.ents.is_empty() {
        return Err("no STEP entities found — is this an ISO-10303-21 file?".into());
    }
    let mut b = Builder::new(&doc, opts);
    b.resolve_units();
    b.resolve_geometry();
    b.build_topology();
    b.mesh_all();
    Ok(b.finish())
}

// ────────────────────────────── STEP text → entities ───────────────────────

/// One instance line. A *complex* instance — `#3 = ( A(..) B(..) )`, which is
/// how rational B-splines and unit contexts are written — keeps its parts in
/// `sub` and has an empty `name`.
#[derive(Debug)]
struct Ent {
    name: String,
    args: Vec<String>,
    sub: Vec<(String, Vec<String>)>,
}

struct Doc {
    ents: HashMap<u32, Ent>,
}

impl Doc {
    fn parse(text: &str) -> Self {
        let clean = strip_comments(text);
        let mut ents = HashMap::new();

        // Statements end at a ';' that is outside quotes and outside parens.
        let mut stmt = String::new();
        let mut depth = 0i32;
        let mut quoted = false;
        for c in clean.chars() {
            match c {
                '\'' => {
                    quoted = !quoted;
                    stmt.push(c);
                }
                '(' if !quoted => {
                    depth += 1;
                    stmt.push(c);
                }
                ')' if !quoted => {
                    depth -= 1;
                    stmt.push(c);
                }
                ';' if !quoted && depth <= 0 => {
                    Self::take_statement(&mut ents, &stmt);
                    stmt.clear();
                }
                _ => stmt.push(c),
            }
        }
        Self { ents }
    }

    fn take_statement(ents: &mut HashMap<u32, Ent>, stmt: &str) {
        let s = stmt.trim();
        if !s.starts_with('#') {
            return;
        }
        let eq = match s.find('=') {
            Some(i) => i,
            None => return,
        };
        let id: u32 = match s[1..eq].trim().parse() {
            Ok(v) => v,
            Err(_) => return,
        };
        let rhs = s[eq + 1..].trim();

        if rhs.starts_with('(') {
            // Complex instance: a parenthesised run of NAME(args) items.
            let inner = &rhs[1..rhs.len().saturating_sub(1)];
            let mut sub = Vec::new();
            for part in split_args_generic(inner, ' ') {
                if let Some(p) = split_call(&part) {
                    sub.push(p);
                }
            }
            ents.insert(
                id,
                Ent {
                    name: String::new(),
                    args: Vec::new(),
                    sub,
                },
            );
        } else if let Some((name, args)) = split_call(rhs) {
            ents.insert(
                id,
                Ent {
                    name,
                    args,
                    sub: Vec::new(),
                },
            );
        }
    }

    /// Type name of an instance; for a complex instance the first part that is
    /// not a bare marker like `GEOMETRIC_REPRESENTATION_ITEM()`.
    fn type_of(&self, id: u32) -> &str {
        match self.ents.get(&id) {
            None => "",
            Some(e) if !e.name.is_empty() => &e.name,
            Some(e) => e
                .sub
                .iter()
                .map(|(n, _)| n.as_str())
                .find(|n| !matches!(*n, "REPRESENTATION_ITEM" | "GEOMETRIC_REPRESENTATION_ITEM" | "CURVE" | "SURFACE" | "BOUNDED_CURVE" | "BOUNDED_SURFACE"))
                .unwrap_or(""),
        }
    }

    /// Arguments of a named entity, looking inside complex instances.
    fn part(&self, id: u32, name: &str) -> Option<&[String]> {
        let e = self.ents.get(&id)?;
        if e.name == name {
            return Some(&e.args);
        }
        e.sub
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, a)| a.as_slice())
    }

    fn has(&self, id: u32, name: &str) -> bool {
        self.part(id, name).is_some()
    }

    fn ids_of(&self, name: &str) -> Vec<u32> {
        let mut v: Vec<u32> = self
            .ents
            .iter()
            .filter(|(_, e)| e.name == name)
            .map(|(&i, _)| i)
            .collect();
        v.sort_unstable();
        v
    }
}

fn strip_comments(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let bytes: Vec<char> = text.chars().collect();
    let mut i = 0;
    let mut quoted = false;
    while i < bytes.len() {
        let c = bytes[i];
        if quoted {
            out.push(c);
            if c == '\'' {
                quoted = false;
            }
            i += 1;
        } else if c == '\'' {
            quoted = true;
            out.push(c);
            i += 1;
        } else if c == '/' && i + 1 < bytes.len() && bytes[i + 1] == '*' {
            i += 2;
            while i + 1 < bytes.len() && !(bytes[i] == '*' && bytes[i + 1] == '/') {
                i += 1;
            }
            i += 2;
        } else {
            out.push(c);
            i += 1;
        }
    }
    out
}

/// Split a `NAME(args)` call into its two halves.
fn split_call(s: &str) -> Option<(String, Vec<String>)> {
    let s = s.trim();
    let open = s.find('(')?;
    if !s.ends_with(')') {
        return None;
    }
    let name = s[..open].trim().to_string();
    if name.is_empty() {
        return None;
    }
    Some((name, split_args(&s[open + 1..s.len() - 1])))
}

/// Split a comma-separated argument list, respecting nesting and quotes.
fn split_args(s: &str) -> Vec<String> {
    split_args_generic(s, ',')
}

/// As [`split_args`], but splitting on `sep` at depth zero. A space separator
/// collapses runs, which is what a complex-instance body needs.
fn split_args_generic(s: &str, sep: char) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut depth = 0i32;
    let mut quoted = false;
    for c in s.chars() {
        match c {
            '\'' => {
                quoted = !quoted;
                cur.push(c);
            }
            '(' if !quoted => {
                depth += 1;
                cur.push(c);
            }
            ')' if !quoted => {
                depth -= 1;
                cur.push(c);
                // A space-separated complex body closes an item at depth 0.
                if sep == ' ' && depth == 0 {
                    out.push(cur.trim().to_string());
                    cur.clear();
                }
            }
            c if c == sep && !quoted && depth == 0 => {
                if sep != ' ' {
                    out.push(cur.trim().to_string());
                    cur.clear();
                }
            }
            _ => cur.push(c),
        }
    }
    if !cur.trim().is_empty() && sep != ' ' {
        out.push(cur.trim().to_string());
    }
    out
}

fn as_ref(s: &str) -> Option<u32> {
    let t = s.trim();
    t.strip_prefix('#')?.trim().parse().ok()
}

fn as_ref_list(s: &str) -> Vec<u32> {
    let t = s.trim();
    let inner = t.strip_prefix('(').and_then(|x| x.strip_suffix(')')).unwrap_or(t);
    split_args(inner).iter().filter_map(|x| as_ref(x)).collect()
}

fn as_f64(s: &str) -> Option<f64> {
    s.trim().parse::<f64>().ok()
}

fn as_f64_list(s: &str) -> Vec<f64> {
    let t = s.trim();
    let inner = t.strip_prefix('(').and_then(|x| x.strip_suffix(')')).unwrap_or(t);
    split_args(inner).iter().filter_map(|x| as_f64(x)).collect()
}

fn as_usize_list(s: &str) -> Vec<usize> {
    let t = s.trim();
    let inner = t.strip_prefix('(').and_then(|x| x.strip_suffix(')')).unwrap_or(t);
    split_args(inner)
        .iter()
        .filter_map(|x| x.trim().parse::<f64>().ok())
        .map(|v| v.max(0.0) as usize)
        .collect()
}

fn as_bool(s: &str) -> bool {
    s.trim() == ".T."
}

// ────────────────────────────────── curves ─────────────────────────────────

#[derive(Clone, Debug)]
struct BSpline {
    degree: usize,
    ctrl: Vec<Point3>,
    weights: Vec<f64>,
    /// Flattened knot vector (multiplicities already expanded).
    knots: Vec<f64>,
}

impl BSpline {
    fn domain(&self) -> (f64, f64) {
        let d = self.degree;
        (self.knots[d], self.knots[self.knots.len() - 1 - d])
    }

    fn eval(&self, t: f64) -> Point3 {
        let h = self.eval_h(t);
        let w = if h[3].abs() < 1e-12 { 1.0 } else { h[3] };
        Point3::new(h[0] / w, h[1] / w, h[2] / w)
    }

    /// de Boor evaluation in homogeneous coordinates, so rational and
    /// non-rational curves take the same path — and so a surface can blend
    /// rows without dividing through twice, which would not be the same curve.
    fn eval_h(&self, t: f64) -> [f64; 4] {
        let (lo, hi) = self.domain();
        let t = t.clamp(lo, hi);
        let p = self.degree;
        let n = self.ctrl.len();

        // Knot span: the last k with knots[k] <= t < knots[k+1], clamped so the
        // final point evaluates on the last span rather than off the end.
        let mut k = p;
        while k + 1 < n && self.knots[k + 1] <= t {
            k += 1;
        }

        let mut dx = Vec::with_capacity(p + 1);
        for j in 0..=p {
            let i = (k + j).saturating_sub(p).min(n - 1);
            let w = self.weights[i];
            dx.push([
                self.ctrl[i].x * w,
                self.ctrl[i].y * w,
                self.ctrl[i].z * w,
                w,
            ]);
        }
        for r in 1..=p {
            for j in (r..=p).rev() {
                let i = k + j - p;
                let lo_k = self.knots[i];
                let hi_k = self.knots.get(i + p + 1 - r).copied().unwrap_or(lo_k + 1.0);
                let den = hi_k - lo_k;
                let a = if den.abs() < 1e-12 { 0.0 } else { (t - lo_k) / den };
                for c in 0..4 {
                    dx[j][c] = (1.0 - a) * dx[j - 1][c] + a * dx[j][c];
                }
            }
        }
        dx[p]
    }
}

#[derive(Clone, Debug)]
enum Curve {
    Line { origin: Point3, dir: Vec3 },
    Circle { frame: Frame3, r: f64 },
    Ellipse { frame: Frame3, a: f64, b: f64 },
    Poly(Vec<Point3>),
    Spline(BSpline),
}

impl Curve {
    fn eval(&self, t: f64) -> Point3 {
        match self {
            Curve::Line { origin, dir } => *origin + *dir * t,
            Curve::Circle { frame, r } => {
                let (s, c) = t.sin_cos();
                frame.origin + frame.x * (c * r) + frame.y * (s * r)
            }
            Curve::Ellipse { frame, a, b } => {
                let (s, c) = t.sin_cos();
                frame.origin + frame.x * (c * a) + frame.y * (s * b)
            }
            Curve::Poly(pts) => {
                if pts.is_empty() {
                    return Point3::ORIGIN;
                }
                let last = (pts.len() - 1) as f64;
                let t = t.clamp(0.0, last);
                let i = t.floor() as usize;
                let f = t - i as f64;
                if i + 1 >= pts.len() {
                    pts[pts.len() - 1]
                } else {
                    pts[i].lerp(pts[i + 1], f)
                }
            }
            Curve::Spline(s) => s.eval(t),
        }
    }

    /// The period, for curves that close on themselves. Only these need the
    /// edge's sense flag to decide *which* arc between two vertices is meant.
    fn period(&self) -> Option<f64> {
        match self {
            Curve::Circle { .. } | Curve::Ellipse { .. } => Some(TAU),
            _ => None,
        }
    }

    /// Parameter of the point on the curve nearest `p`.
    fn param_of(&self, p: Point3) -> f64 {
        match self {
            Curve::Line { origin, dir } => {
                let l2 = dir.length_sq();
                if l2 < 1e-18 {
                    0.0
                } else {
                    (p - *origin).dot(*dir) / l2
                }
            }
            Curve::Circle { frame, .. } => {
                let l = frame.to_local_point(p);
                l.y.atan2(l.x)
            }
            Curve::Ellipse { frame, a, b } => {
                let l = frame.to_local_point(p);
                (l.y / b.max(1e-12)).atan2(l.x / a.max(1e-12))
            }
            Curve::Poly(pts) => {
                let mut best = 0.0;
                let mut best_d = f64::INFINITY;
                for (i, q) in pts.iter().enumerate() {
                    let d = (p - *q).length_sq();
                    if d < best_d {
                        best_d = d;
                        best = i as f64;
                    }
                }
                best
            }
            Curve::Spline(s) => {
                // Coarse scan then a short golden-section polish: the vertices
                // of an imported edge sit on the curve, so this converges fast
                // and only ever has to be right to within the chord tolerance.
                let (lo, hi) = s.domain();
                let n = 256;
                let mut best = lo;
                let mut best_d = f64::INFINITY;
                for i in 0..=n {
                    let t = lo + (hi - lo) * i as f64 / n as f64;
                    let d = (s.eval(t) - p).length_sq();
                    if d < best_d {
                        best_d = d;
                        best = t;
                    }
                }
                let step = (hi - lo) / n as f64;
                let (mut a, mut b) = ((best - step).max(lo), (best + step).min(hi));
                for _ in 0..40 {
                    let m1 = a + (b - a) / 3.0;
                    let m2 = b - (b - a) / 3.0;
                    if (s.eval(m1) - p).length_sq() < (s.eval(m2) - p).length_sq() {
                        b = m2;
                    } else {
                        a = m1;
                    }
                }
                0.5 * (a + b)
            }
        }
    }

}

// ───────────────────────────────── surfaces ────────────────────────────────

#[derive(Clone, Debug)]
struct BSplineSurf {
    du: usize,
    dv: usize,
    /// Control net, `ctrl[i][j]` with `i` along u.
    ctrl: Vec<Vec<Point3>>,
    weights: Vec<Vec<f64>>,
    ku: Vec<f64>,
    kv: Vec<f64>,
    /// Closure flags from the file. A NURBS that closes on itself is periodic
    /// in exactly the way a cylinder is, and its faces are bands; missing that
    /// sends every elbow patch down the wrong path.
    closed_u: bool,
    closed_v: bool,
    /// Millimetres per unit of parameter, measured once from the surface.
    mu: f64,
    mv: f64,
    /// A coarse table of surface points, built once, used to seed the
    /// point-to-parameter search. Without it every boundary point of every
    /// patch re-evaluates the whole surface a few hundred times, which is most
    /// of the time a large model takes to open.
    seed: Vec<(f64, f64, Point3)>,
}

impl BSplineSurf {
    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        (
            (self.ku[self.du], self.ku[self.ku.len() - 1 - self.du]),
            (self.kv[self.dv], self.kv[self.kv.len() - 1 - self.dv]),
        )
    }

    /// Evaluate the v-isocurve of every u-row in homogeneous coordinates, then
    /// blend those along u with the weights those rows came out with. Carrying
    /// the weight through is what makes a rational patch — a quarter-torus
    /// written as a NURBS, say — come back as the exact quadric it encodes.
    fn eval(&self, u: f64, v: f64) -> Point3 {
        // A closed patch is periodic, and the mesher works in unwrapped
        // parameters: a rim that crosses the seam carries on past the end of
        // the knot vector (v = -2.7 on a tube parametrised 0..4). Outside its
        // knots a B-spline just clamps to its end point, so without folding the
        // parameter back every interior row of such a band collapsed onto one
        // point, and each U-turn of a serpentine came back as a flat chord
        // across its bend with a fifth of the tube gone.
        let ((u0, u1), (v0, v1)) = self.domain();
        let u = if self.closed_u && u1 > u0 { u0 + (u - u0).rem_euclid(u1 - u0) } else { u };
        let v = if self.closed_v && v1 > v0 { v0 + (v - v0).rem_euclid(v1 - v0) } else { v };
        let n = self.ctrl.len();
        let mut pts = Vec::with_capacity(n);
        let mut ws = Vec::with_capacity(n);
        for i in 0..n {
            let c = BSpline {
                degree: self.dv,
                ctrl: self.ctrl[i].clone(),
                weights: self.weights[i].clone(),
                knots: self.kv.clone(),
            };
            let h = c.eval_h(v);
            let w = if h[3].abs() < 1e-12 { 1.0 } else { h[3] };
            pts.push(Point3::new(h[0] / w, h[1] / w, h[2] / w));
            ws.push(w);
        }
        BSpline {
            degree: self.du,
            ctrl: pts,
            weights: ws,
            knots: self.ku.clone(),
        }
        .eval(u)
    }
}

#[derive(Clone, Debug)]
enum Surf {
    Plane(Frame3),
    Cyl { frame: Frame3, r: f64 },
    Cone { frame: Frame3, r: f64, half_angle: f64 },
    Sphere { frame: Frame3, r: f64 },
    Torus { frame: Frame3, major: f64, minor: f64 },
    Spline(BSplineSurf),
}

impl Surf {
    fn eval(&self, u: f64, v: f64) -> Point3 {
        match self {
            Surf::Plane(f) => f.origin + f.x * u + f.y * v,
            Surf::Cyl { frame, r } => {
                let (s, c) = u.sin_cos();
                frame.origin + frame.x * (c * r) + frame.y * (s * r) + frame.z * v
            }
            Surf::Cone { frame, r, half_angle } => {
                let rr = r + v * half_angle.tan();
                let (s, c) = u.sin_cos();
                frame.origin + frame.x * (c * rr) + frame.y * (s * rr) + frame.z * v
            }
            Surf::Sphere { frame, r } => {
                let (su, cu) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                frame.origin + frame.x * (cv * cu * r) + frame.y * (cv * su * r) + frame.z * (sv * r)
            }
            Surf::Torus { frame, major, minor } => {
                let (su, cu) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                let rr = major + minor * cv;
                frame.origin + frame.x * (cu * rr) + frame.y * (su * rr) + frame.z * (minor * sv)
            }
            Surf::Spline(s) => s.eval(u, v),
        }
    }

    /// The natural normal `∂P/∂u × ∂P/∂v`, normalised. For every analytic case
    /// here that direction points out of the material, which is what STEP's
    /// `same_sense` flag is stated relative to.
    fn normal(&self, u: f64, v: f64) -> Vec3 {
        match self {
            Surf::Plane(f) => f.z.as_vec(),
            Surf::Cyl { frame, .. } => {
                let (s, c) = u.sin_cos();
                frame.x * c + frame.y * s
            }
            Surf::Sphere { frame, .. } => {
                let (su, cu) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                frame.x * (cv * cu) + frame.y * (cv * su) + frame.z * sv
            }
            Surf::Torus { frame, .. } => {
                let (su, cu) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                (frame.x * cu + frame.y * su) * cv + frame.z.as_vec() * sv
            }
            // Cone and spline: finite differences, guaranteed consistent with
            // `eval` even where the closed form would be fiddly. The step is
            // taken in parameter units that amount to a fixed distance on the
            // surface, so a patch parametrised 0..1 and one parametrised 0..2pi
            // get the same answer.
            _ => {
                let (mu, mv) = self.metric();
                let (hu, hv) = (1e-4 / mu, 1e-4 / mv);
                let du = self.eval(u + hu, v) - self.eval(u - hu, v);
                let dv = self.eval(u, v + hv) - self.eval(u, v - hv);
                match (du.try_normalize(), dv.try_normalize()) {
                    (Some(a), Some(b)) => a.cross(b).try_normalize().unwrap_or(Vec3::Z),
                    _ => Vec3::Z,
                }
            }
        }
    }

    fn uv_of(&self, p: Point3) -> (f64, f64) {
        match self {
            Surf::Plane(f) => {
                let l = f.to_local_point(p);
                (l.x, l.y)
            }
            Surf::Cyl { frame, .. } | Surf::Cone { frame, .. } => {
                let l = frame.to_local_point(p);
                (l.y.atan2(l.x), l.z)
            }
            Surf::Sphere { frame, r } => {
                let l = frame.to_local_point(p);
                ((l.y).atan2(l.x), (l.z / r.max(1e-12)).clamp(-1.0, 1.0).asin())
            }
            Surf::Torus { frame, major, .. } => {
                let l = frame.to_local_point(p);
                let theta = l.y.atan2(l.x);
                let radial = (l.x * l.x + l.y * l.y).sqrt() - major;
                (theta, l.z.atan2(radial))
            }
            Surf::Spline(s) => {
                let ((u0, u1), (v0, v1)) = s.domain();
                let (mut u, mut v) = (u0, v0);
                let mut best = f64::INFINITY;
                for &(su, sv, sp) in &s.seed {
                    let d = (sp - p).length_sq();
                    if d < best {
                        best = d;
                        u = su;
                        v = sv;
                    }
                }
                // Local descent from the nearest table entry; a handful of
                // halvings is plenty when the seed is already within one cell.
                let n = (s.seed.len() as f64).sqrt().max(2.0) - 1.0;
                let mut hu = (u1 - u0) / n;
                let mut hv = (v1 - v0) / n;
                for _ in 0..30 {
                    let mut improved = false;
                    for (du, dv) in [(hu, 0.0), (-hu, 0.0), (0.0, hv), (0.0, -hv)] {
                        let (cu, cv) = ((u + du).clamp(u0, u1), (v + dv).clamp(v0, v1));
                        let d = (s.eval(cu, cv) - p).length_sq();
                        if d < best {
                            best = d;
                            u = cu;
                            v = cv;
                            improved = true;
                        }
                    }
                    if !improved {
                        hu *= 0.5;
                        hv *= 0.5;
                    }
                }
                (u, v)
            }
        }
    }

    /// Period of each parameter, where it wraps.
    fn periods(&self) -> (Option<f64>, Option<f64>) {
        match self {
            Surf::Plane(_) => (None, None),
            Surf::Cyl { .. } | Surf::Cone { .. } | Surf::Sphere { .. } => (Some(TAU), None),
            Surf::Torus { .. } => (Some(TAU), Some(TAU)),
            Surf::Spline(s) => {
                let ((u0, u1), (v0, v1)) = s.domain();
                (
                    s.closed_u.then_some(u1 - u0),
                    s.closed_v.then_some(v1 - v0),
                )
            }
        }
    }

    /// Scale factors turning parameter steps into roughly millimetres, so that
    /// triangulation and refinement behave the same on a 0.3 mm filament as on
    /// a 15 mm plate.
    fn metric(&self) -> (f64, f64) {
        match self {
            Surf::Plane(_) => (1.0, 1.0),
            Surf::Cyl { r, .. } => (r.max(1e-6), 1.0),
            Surf::Cone { r, .. } => (r.abs().max(1e-6), 1.0),
            Surf::Sphere { r, .. } => (r.max(1e-6), r.max(1e-6)),
            Surf::Torus { major, minor, .. } => (major.max(1e-6), minor.max(1e-6)),
            Surf::Spline(s) => (s.mu, s.mv),
        }
    }

    /// True when straight lines in `(u, v)` leave the surface, i.e. when
    /// interior refinement buys anything.
    fn is_curved(&self) -> bool {
        !matches!(self, Surf::Plane(_))
    }

    /// `v` values at which the surface collapses to a point, so a boundary that
    /// wraps once can be closed against the pole.
    fn poles(&self) -> Vec<f64> {
        match self {
            Surf::Sphere { .. } => vec![-PI / 2.0, PI / 2.0],
            Surf::Cone { r, half_angle, .. } => {
                let t = half_angle.tan();
                if t.abs() < 1e-9 {
                    vec![]
                } else {
                    vec![-r / t]
                }
            }
            _ => vec![],
        }
    }
}

// ─────────────────────────────── the builder ───────────────────────────────

struct FaceIn {
    surf_id: u32,
    same_sense: bool,
    /// Loops as lists of `(edge id, forward)`; a vertex loop is empty.
    loops: Vec<Vec<(u32, bool)>>,
    solid: u32,
}

struct Builder<'a> {
    doc: &'a Doc,
    opts: ImportOptions,
    scale: f64,

    points: HashMap<u32, Point3>,
    dirs: HashMap<u32, Vec3>,
    frames: HashMap<u32, Frame3>,
    curves: HashMap<u32, Curve>,
    surfs: HashMap<u32, Surf>,

    /// `EDGE_CURVE` id → (start vertex point, end vertex point, curve id, sense)
    edges: HashMap<u32, (Point3, Point3, u32, bool)>,
    /// Discretised edges, cached so both faces of an edge get identical points.
    edge_pts: HashMap<u32, Vec<Point3>>,

    faces: Vec<(u32, FaceIn)>,

    out: Vec<FaceMesh>,
    report: ImportReport,
    fail_reasons: HashMap<String, usize>,
    warn_set: HashMap<String, usize>,
}

impl<'a> Builder<'a> {
    fn new(doc: &'a Doc, opts: &ImportOptions) -> Self {
        Self {
            doc,
            opts: opts.clone(),
            scale: 1.0,
            points: HashMap::new(),
            dirs: HashMap::new(),
            frames: HashMap::new(),
            curves: HashMap::new(),
            surfs: HashMap::new(),
            edges: HashMap::new(),
            edge_pts: HashMap::new(),
            faces: Vec::new(),
            out: Vec::new(),
            report: ImportReport::default(),
            fail_reasons: HashMap::new(),
            warn_set: HashMap::new(),
        }
    }

    fn warn(&mut self, msg: impl Into<String>) {
        *self.warn_set.entry(msg.into()).or_insert(0) += 1;
    }

    fn fail(&mut self, msg: impl Into<String>) {
        *self.fail_reasons.entry(msg.into()).or_insert(0) += 1;
        self.report.faces_failed += 1;
    }

    // ── units ──────────────────────────────────────────────────────────────

    /// Millimetres per file unit, taken from the file's own unit context.
    ///
    /// Getting this wrong scales the whole model, so it is read rather than
    /// assumed; when the file does not say, that fact is reported instead of
    /// being hidden behind a default.
    fn resolve_units(&mut self) {
        let mut found: Option<f64> = None;
        for (&id, e) in &self.doc.ents {
            let is_len = e.name == "SI_UNIT"
                || e.sub.iter().any(|(n, _)| n == "LENGTH_UNIT")
                || self.doc.has(id, "LENGTH_UNIT");
            if !is_len {
                continue;
            }
            if let Some(args) = self.doc.part(id, "SI_UNIT") {
                if args.len() >= 2 && args[1].trim() == ".METRE." {
                    let mm = match args[0].trim() {
                        ".MILLI." => 1.0,
                        ".CENTI." => 10.0,
                        ".DECI." => 100.0,
                        "$" => 1000.0,
                        ".KILO." => 1.0e6,
                        ".MICRO." => 1.0e-3,
                        _ => 1.0,
                    };
                    found = Some(mm);
                    // A millimetre unit is by far the most common; keep looking
                    // only if this was a non-length SI unit.
                    if (mm - 1.0).abs() < 1e-12 {
                        break;
                    }
                }
            }
        }
        // Inch and other conversion-based units.
        if found.is_none() {
            for id in self.doc.ids_of("CONVERSION_BASED_UNIT") {
                if let Some(args) = self.doc.part(id, "CONVERSION_BASED_UNIT") {
                    if args.len() >= 2 {
                        if let Some(mf) = as_ref(&args[1]) {
                            if let Some(m) = self.doc.part(mf, "MEASURE_WITH_UNIT") {
                                if let Some(v) = m.first().and_then(|s| {
                                    s.split('(').nth(1).and_then(|x| as_f64(x.trim_end_matches(')')))
                                }) {
                                    found = Some(v * 1000.0);
                                }
                            }
                        }
                    }
                }
            }
        }
        match found {
            Some(mm) => self.scale = mm,
            None => {
                self.scale = 1.0;
                self.warn("file declares no length unit; assuming millimetres");
            }
        }
        self.report.length_scale_mm = self.scale;
    }

    // ── geometry ───────────────────────────────────────────────────────────

    fn resolve_geometry(&mut self) {
        let s = self.scale;
        for (&id, e) in &self.doc.ents {
            match e.name.as_str() {
                "CARTESIAN_POINT" => {
                    if let Some(a) = e.args.get(1) {
                        let c = as_f64_list(a);
                        if c.len() >= 3 {
                            self.points.insert(id, Point3::new(c[0] * s, c[1] * s, c[2] * s));
                        }
                    }
                }
                "DIRECTION" => {
                    if let Some(a) = e.args.get(1) {
                        let c = as_f64_list(a);
                        if c.len() >= 3 {
                            self.dirs.insert(id, Vec3::new(c[0], c[1], c[2]));
                        }
                    }
                }
                _ => {}
            }
        }

        // VECTOR wraps a direction with a magnitude; LINE points at one of
        // these, not at a bare DIRECTION. The old reader missed that and turned
        // every straight edge into an X-axis line.
        let mut vectors: HashMap<u32, Vec3> = HashMap::new();
        for id in self.doc.ids_of("VECTOR") {
            if let Some(a) = self.doc.part(id, "VECTOR") {
                if a.len() >= 3 {
                    let d = as_ref(&a[1]).and_then(|r| self.dirs.get(&r).copied()).unwrap_or(Vec3::X);
                    let m = as_f64(&a[2]).unwrap_or(1.0) * s;
                    if let Some(u) = d.try_normalize() {
                        vectors.insert(id, u * m);
                    }
                }
            }
        }

        for id in self.doc.ids_of("AXIS2_PLACEMENT_3D") {
            if let Some(f) = self.axis2(id) {
                self.frames.insert(id, f);
            }
        }

        // Curves. Dispatch is by *which part an instance carries*, not by its
        // type name: a rational B-spline arrives as a complex instance whose
        // leading part is `B_SPLINE_CURVE`, so matching on the name alone
        // silently drops every NURBS in the file.
        let ids: Vec<u32> = self.doc.ents.keys().copied().collect();
        for id in ids {
            let c = if self.doc.has(id, "B_SPLINE_CURVE_WITH_KNOTS") {
                self.bspline_curve(id).map(Curve::Spline)
            } else {
                let ty = self.doc.type_of(id).to_string();
                match ty.as_str() {
                "LINE" => self.doc.part(id, "LINE").and_then(|a| {
                    let o = as_ref(a.get(1)?).and_then(|r| self.points.get(&r).copied())?;
                    let d = as_ref(a.get(2)?)
                        .and_then(|r| vectors.get(&r).copied().or_else(|| self.dirs.get(&r).copied()))?;
                    Some(Curve::Line { origin: o, dir: d })
                }),
                "CIRCLE" => self.doc.part(id, "CIRCLE").and_then(|a| {
                    let f = as_ref(a.get(1)?).and_then(|r| self.frames.get(&r).copied())?;
                    Some(Curve::Circle {
                        frame: f,
                        r: as_f64(a.get(2)?)? * s,
                    })
                }),
                "ELLIPSE" => self.doc.part(id, "ELLIPSE").and_then(|a| {
                    let f = as_ref(a.get(1)?).and_then(|r| self.frames.get(&r).copied())?;
                    Some(Curve::Ellipse {
                        frame: f,
                        a: as_f64(a.get(2)?)? * s,
                        b: as_f64(a.get(3)?)? * s,
                    })
                }),
                "POLYLINE" => self.doc.part(id, "POLYLINE").and_then(|a| {
                    let pts: Vec<Point3> = as_ref_list(a.get(1)?)
                        .iter()
                        .filter_map(|r| self.points.get(r).copied())
                        .collect();
                    (pts.len() >= 2).then_some(Curve::Poly(pts))
                }),
                _ => None,
                }
            };
            if let Some(c) = c {
                self.curves.insert(id, c);
            }
        }

        // Surfaces, same rule.
        let ids: Vec<u32> = self.doc.ents.keys().copied().collect();
        for id in ids {
            let sf = if self.doc.has(id, "B_SPLINE_SURFACE_WITH_KNOTS") {
                self.bspline_surface(id).map(Surf::Spline)
            } else {
                let ty = self.doc.type_of(id).to_string();
                match ty.as_str() {
                "PLANE" => self
                    .doc
                    .part(id, "PLANE")
                    .and_then(|a| as_ref(a.get(1)?))
                    .and_then(|r| self.frames.get(&r).copied())
                    .map(Surf::Plane),
                "CYLINDRICAL_SURFACE" => self.doc.part(id, "CYLINDRICAL_SURFACE").and_then(|a| {
                    let f = as_ref(a.get(1)?).and_then(|r| self.frames.get(&r).copied())?;
                    Some(Surf::Cyl {
                        frame: f,
                        r: as_f64(a.get(2)?)? * s,
                    })
                }),
                "CONICAL_SURFACE" => self.doc.part(id, "CONICAL_SURFACE").and_then(|a| {
                    let f = as_ref(a.get(1)?).and_then(|r| self.frames.get(&r).copied())?;
                    Some(Surf::Cone {
                        frame: f,
                        r: as_f64(a.get(2)?)? * s,
                        half_angle: as_f64(a.get(3)?)?,
                    })
                }),
                "SPHERICAL_SURFACE" => self.doc.part(id, "SPHERICAL_SURFACE").and_then(|a| {
                    let f = as_ref(a.get(1)?).and_then(|r| self.frames.get(&r).copied())?;
                    Some(Surf::Sphere {
                        frame: f,
                        r: as_f64(a.get(2)?)? * s,
                    })
                }),
                "TOROIDAL_SURFACE" => self.doc.part(id, "TOROIDAL_SURFACE").and_then(|a| {
                    let f = as_ref(a.get(1)?).and_then(|r| self.frames.get(&r).copied())?;
                    Some(Surf::Torus {
                        frame: f,
                        major: as_f64(a.get(2)?)? * s,
                        minor: as_f64(a.get(3)?)? * s,
                    })
                }),
                _ => None,
                }
            };
            if let Some(sf) = sf {
                self.surfs.insert(id, sf);
            }
        }
    }

    /// `AXIS2_PLACEMENT_3D` with the reference direction properly projected off
    /// the axis. STEP does not require the two to be perpendicular, and a
    /// non-orthogonal frame rotates the surface parametrisation, which then
    /// trims the wrong arc.
    fn axis2(&self, id: u32) -> Option<Frame3> {
        let a = self.doc.part(id, "AXIS2_PLACEMENT_3D")?;
        let origin = as_ref(a.get(1)?)
            .and_then(|r| self.points.get(&r).copied())
            .unwrap_or(Point3::ORIGIN);
        let z = a
            .get(2)
            .and_then(|s| as_ref(s))
            .and_then(|r| self.dirs.get(&r).copied())
            .and_then(UnitVec3::try_from_vec)
            .unwrap_or(UnitVec3::Z);
        let x_ref = a.get(3).and_then(|s| as_ref(s)).and_then(|r| self.dirs.get(&r).copied());
        let x = match x_ref {
            Some(v) => {
                let proj = v - z.as_vec() * z.dot_vec(v);
                UnitVec3::try_from_vec(proj).unwrap_or_else(|| z.perp_basis().0)
            }
            None => z.perp_basis().0,
        };
        let y = UnitVec3::try_from_vec(z.cross(x)).unwrap_or_else(|| z.perp_basis().1);
        Some(Frame3 { origin, x, y, z })
    }

    fn bspline_curve(&self, id: u32) -> Option<BSpline> {
        let base = self.doc.part(id, "B_SPLINE_CURVE").or_else(|| self.doc.part(id, "B_SPLINE_CURVE_WITH_KNOTS"))?;
        let degree = as_f64(base.first()?)? as usize;
        let ctrl: Vec<Point3> = as_ref_list(base.get(1)?)
            .iter()
            .filter_map(|r| self.points.get(r).copied())
            .collect();
        if ctrl.len() < degree + 1 {
            return None;
        }
        // In a complex instance the knots live in the WITH_KNOTS part, which
        // then carries only the knot data.
        let kargs = self.doc.part(id, "B_SPLINE_CURVE_WITH_KNOTS")?;
        let (mult_s, knot_s) = if kargs.len() >= 6 {
            (kargs.get(5)?, kargs.get(6)?)
        } else {
            (kargs.first()?, kargs.get(1)?)
        };
        let mults = as_usize_list(mult_s);
        let kn = as_f64_list(knot_s);
        if mults.len() != kn.len() || kn.is_empty() {
            return None;
        }
        let mut knots = Vec::new();
        for (m, k) in mults.iter().zip(kn.iter()) {
            for _ in 0..*m {
                knots.push(*k);
            }
        }
        let weights = self
            .doc
            .part(id, "RATIONAL_B_SPLINE_CURVE")
            .and_then(|a| a.first().map(|s| as_f64_list(s)))
            .filter(|w| w.len() == ctrl.len())
            .unwrap_or_else(|| vec![1.0; ctrl.len()]);
        Some(BSpline {
            degree,
            ctrl,
            weights,
            knots,
        })
    }

    fn bspline_surface(&self, id: u32) -> Option<BSplineSurf> {
        let base = self
            .doc
            .part(id, "B_SPLINE_SURFACE")
            .or_else(|| self.doc.part(id, "B_SPLINE_SURFACE_WITH_KNOTS"))?;
        let du = as_f64(base.first()?)? as usize;
        let dv = as_f64(base.get(1)?)? as usize;
        // Control net: a list of lists of point refs.
        let net_s = base.get(2)?.trim();
        let rows = split_args(net_s.strip_prefix('(')?.strip_suffix(')')?);
        let mut ctrl = Vec::new();
        for r in &rows {
            let pts: Vec<Point3> = as_ref_list(r)
                .iter()
                .filter_map(|x| self.points.get(x).copied())
                .collect();
            if pts.is_empty() {
                return None;
            }
            ctrl.push(pts);
        }
        if ctrl.len() < du + 1 || ctrl[0].len() < dv + 1 {
            return None;
        }
        let ka = self.doc.part(id, "B_SPLINE_SURFACE_WITH_KNOTS")?;
        // Simple form: args 0..6 are degrees, net, form, flags; knots start at 7.
        let off = if ka.len() >= 12 { 7 } else { 0 };
        let mu = as_usize_list(ka.get(off)?);
        let mv = as_usize_list(ka.get(off + 1)?);
        let ku_v = as_f64_list(ka.get(off + 2)?);
        let kv_v = as_f64_list(ka.get(off + 3)?);
        if mu.len() != ku_v.len() || mv.len() != kv_v.len() {
            return None;
        }
        let expand = |m: &[usize], k: &[f64]| {
            let mut o = Vec::new();
            for (mm, kk) in m.iter().zip(k.iter()) {
                for _ in 0..*mm {
                    o.push(*kk);
                }
            }
            o
        };
        let weights = self
            .doc
            .part(id, "RATIONAL_B_SPLINE_SURFACE")
            .and_then(|a| {
                let rows = split_args(a.first()?.trim().strip_prefix('(')?.strip_suffix(')')?);
                let w: Vec<Vec<f64>> = rows.iter().map(|r| as_f64_list(r)).collect();
                (w.len() == ctrl.len() && w.iter().zip(&ctrl).all(|(a, b)| a.len() == b.len())).then_some(w)
            })
            .unwrap_or_else(|| ctrl.iter().map(|r| vec![1.0; r.len()]).collect());

        let closed_u = base.get(4).map(|s| as_bool(s)).unwrap_or(false);
        let closed_v = base.get(5).map(|s| as_bool(s)).unwrap_or(false);

        let mut sf = BSplineSurf {
            du,
            dv,
            ctrl,
            weights,
            ku: expand(&mu, &ku_v),
            kv: expand(&mv, &kv_v),
            closed_u,
            closed_v,
            mu: 1.0,
            mv: 1.0,
            seed: Vec::new(),
        };
        // Measure the parametric scale so tolerances mean millimetres here too.
        let ((u0, u1), (v0, v1)) = sf.domain();
        let mid_u = 0.5 * (u0 + u1);
        let mid_v = 0.5 * (v0 + v1);
        let len_u: f64 = (0..8)
            .map(|i| {
                let a = u0 + (u1 - u0) * i as f64 / 8.0;
                let b = u0 + (u1 - u0) * (i + 1) as f64 / 8.0;
                (sf.eval(b, mid_v) - sf.eval(a, mid_v)).length()
            })
            .sum();
        let len_v: f64 = (0..8)
            .map(|i| {
                let a = v0 + (v1 - v0) * i as f64 / 8.0;
                let b = v0 + (v1 - v0) * (i + 1) as f64 / 8.0;
                (sf.eval(mid_u, b) - sf.eval(mid_u, a)).length()
            })
            .sum();
        sf.mu = if u1 > u0 { (len_u / (u1 - u0)).max(1e-6) } else { 1.0 };
        sf.mv = if v1 > v0 { (len_v / (v1 - v0)).max(1e-6) } else { 1.0 };

        const SEED: usize = 20;
        let mut seed = Vec::with_capacity((SEED + 1) * (SEED + 1));
        for i in 0..=SEED {
            for j in 0..=SEED {
                let u = u0 + (u1 - u0) * i as f64 / SEED as f64;
                let v = v0 + (v1 - v0) * j as f64 / SEED as f64;
                seed.push((u, v, sf.eval(u, v)));
            }
        }
        sf.seed = seed;
        Some(sf)
    }

    // ── topology ───────────────────────────────────────────────────────────

    fn build_topology(&mut self) {
        let mut vertices: HashMap<u32, Point3> = HashMap::new();
        for id in self.doc.ids_of("VERTEX_POINT") {
            if let Some(a) = self.doc.part(id, "VERTEX_POINT") {
                if let Some(p) = a.get(1).and_then(|s| as_ref(s)).and_then(|r| self.points.get(&r)) {
                    vertices.insert(id, *p);
                }
            }
        }

        for id in self.doc.ids_of("EDGE_CURVE") {
            if let Some(a) = self.doc.part(id, "EDGE_CURVE") {
                if a.len() < 5 {
                    continue;
                }
                let v1 = as_ref(&a[1]).and_then(|r| vertices.get(&r).copied());
                let v2 = as_ref(&a[2]).and_then(|r| vertices.get(&r).copied());
                let cref = as_ref(&a[3]);
                if let (Some(v1), Some(v2), Some(c)) = (v1, v2, cref) {
                    let c = self.resolve_curve_ref(c);
                    self.edges.insert(id, (v1, v2, c, as_bool(&a[4])));
                }
            }
        }

        // Oriented edges, resolved to (edge, direction).
        let mut oriented: HashMap<u32, (u32, bool)> = HashMap::new();
        for id in self.doc.ids_of("ORIENTED_EDGE") {
            if let Some(a) = self.doc.part(id, "ORIENTED_EDGE") {
                if a.len() >= 5 {
                    if let Some(e) = as_ref(&a[3]) {
                        oriented.insert(id, (e, as_bool(&a[4])));
                    }
                }
            }
        }

        let mut loops: HashMap<u32, Vec<(u32, bool)>> = HashMap::new();
        for id in self.doc.ids_of("EDGE_LOOP") {
            if let Some(a) = self.doc.part(id, "EDGE_LOOP") {
                if let Some(list) = a.get(1) {
                    let mut v = Vec::new();
                    for oe in as_ref_list(list) {
                        if let Some(&(e, fwd)) = oriented.get(&oe) {
                            v.push((e, fwd));
                        }
                    }
                    loops.insert(id, v);
                }
            }
        }
        for id in self.doc.ids_of("VERTEX_LOOP") {
            loops.insert(id, Vec::new());
        }

        // Bounds wrap a loop and may reverse it.
        let mut bounds: HashMap<u32, (u32, bool)> = HashMap::new();
        for name in ["FACE_OUTER_BOUND", "FACE_BOUND"] {
            for id in self.doc.ids_of(name) {
                if let Some(a) = self.doc.part(id, name) {
                    if a.len() >= 3 {
                        if let Some(l) = as_ref(&a[1]) {
                            bounds.insert(id, (l, as_bool(&a[2])));
                        }
                    }
                }
            }
        }

        // Which solid each face belongs to, for the viewer's per-solid grouping.
        let mut face_solid: HashMap<u32, u32> = HashMap::new();
        let mut solids = 0usize;
        for name in ["MANIFOLD_SOLID_BREP", "BREP_WITH_VOIDS", "SHELL_BASED_SURFACE_MODEL"] {
            for id in self.doc.ids_of(name) {
                if let Some(a) = self.doc.part(id, name) {
                    solids += 1;
                    let shells: Vec<u32> = if name == "SHELL_BASED_SURFACE_MODEL" {
                        a.get(1).map(|s| as_ref_list(s)).unwrap_or_default()
                    } else {
                        a.get(1).and_then(|s| as_ref(s)).into_iter().collect()
                    };
                    for sh in shells {
                        let inner = sh;
                        for sname in ["CLOSED_SHELL", "OPEN_SHELL"] {
                            if let Some(sa) = self.doc.part(inner, sname) {
                                if let Some(fl) = sa.get(1) {
                                    for f in as_ref_list(fl) {
                                        face_solid.insert(f, id);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // Shells that no solid claimed still carry faces worth drawing.
        for sname in ["CLOSED_SHELL", "OPEN_SHELL"] {
            for id in self.doc.ids_of(sname) {
                if let Some(sa) = self.doc.part(id, sname) {
                    if let Some(fl) = sa.get(1) {
                        for f in as_ref_list(fl) {
                            face_solid.entry(f).or_insert(id);
                        }
                    }
                }
            }
        }
        self.report.solids = solids.max(1);

        for name in ["ADVANCED_FACE", "FACE_SURFACE"] {
            for id in self.doc.ids_of(name) {
                let a = match self.doc.part(id, name) {
                    Some(a) if a.len() >= 4 => a,
                    _ => continue,
                };
                let bl = as_ref_list(&a[1]);
                let surf_id = match as_ref(&a[2]) {
                    Some(s) => s,
                    None => continue,
                };
                let same_sense = as_bool(&a[3]);
                let mut ls = Vec::new();
                for b in bl {
                    if let Some(&(lp, fwd)) = bounds.get(&b) {
                        if let Some(edges) = loops.get(&lp) {
                            let mut e = edges.clone();
                            if !fwd {
                                e.reverse();
                                for x in &mut e {
                                    x.1 = !x.1;
                                }
                            }
                            ls.push(e);
                        }
                    }
                }
                self.faces.push((
                    id,
                    FaceIn {
                        surf_id,
                        same_sense,
                        loops: ls,
                        solid: face_solid.get(&id).copied().unwrap_or(0),
                    },
                ));
            }
        }
        self.faces.sort_by_key(|(id, _)| *id);
        self.report.faces_total = self.faces.len();

        if !self.doc.ids_of("NEXT_ASSEMBLY_USAGE_OCCURRENCE").is_empty() {
            self.warn("file is an assembly; component placements are not applied");
        }
    }

    /// Follow `SURFACE_CURVE` / `SEAM_CURVE` / `TRIMMED_CURVE` down to the
    /// 3-D curve they are built on.
    fn resolve_curve_ref(&self, mut id: u32) -> u32 {
        for _ in 0..8 {
            let ty = self.doc.type_of(id);
            let next = match ty {
                "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE" | "BOUNDED_SURFACE_CURVE" => self
                    .doc
                    .part(id, ty)
                    .and_then(|a| a.get(1).and_then(|s| as_ref(s))),
                "TRIMMED_CURVE" => self.doc.part(id, ty).and_then(|a| a.first().and_then(|s| as_ref(s))),
                _ => None,
            };
            match next {
                Some(n) if n != id => id = n,
                _ => break,
            }
        }
        id
    }

    // ── meshing ────────────────────────────────────────────────────────────

    fn mesh_all(&mut self) {
        let faces = std::mem::take(&mut self.faces);
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];

        for (fid, f) in &faces {
            let surf = match self.surfs.get(&f.surf_id) {
                Some(s) => s.clone(),
                None => {
                    let ty = self.doc.type_of(f.surf_id).to_string();
                    self.fail(format!(
                        "unsupported surface type {}",
                        if ty.is_empty() { "?" } else { &ty }
                    ));
                    continue;
                }
            };

            // Discretise every loop, reusing cached edge polylines so that the
            // face on the other side of each edge gets the same points.
            let mut rings: Vec<Vec<Point3>> = Vec::new();
            for lp in &f.loops {
                if lp.is_empty() {
                    continue;
                }
                let mut ring: Vec<Point3> = Vec::new();
                for &(e, fwd) in lp {
                    let pts = self.edge_polyline(e);
                    if pts.len() < 2 {
                        continue;
                    }
                    let mut seg: Vec<Point3> = pts.clone();
                    if !fwd {
                        seg.reverse();
                    }
                    // Drop the joint point; the previous edge already placed it.
                    let start = usize::from(!ring.is_empty());
                    ring.extend_from_slice(&seg[start..]);
                }
                // The ring closes back on itself; the duplicate is implicit.
                if ring.len() >= 2 {
                    if let (Some(a), Some(b)) = (ring.first().copied(), ring.last().copied()) {
                        if (a - b).length() < 1e-9 {
                            ring.pop();
                        }
                    }
                }
                if ring.len() >= 3 {
                    rings.push(ring);
                }
            }

            DUMP_FACE.with(|c| c.set(*fid));
            match mesh_face(&surf, &rings, f.same_sense, &self.opts) {
                Ok((pos, nor, idx)) => {
                    if idx.is_empty() {
                        self.fail("face produced no triangles");
                        continue;
                    }

                    for p in &pos {
                        lo[0] = lo[0].min(p.x);
                        lo[1] = lo[1].min(p.y);
                        lo[2] = lo[2].min(p.z);
                        hi[0] = hi[0].max(p.x);
                        hi[1] = hi[1].max(p.y);
                        hi[2] = hi[2].max(p.z);
                    }
                    self.report.faces_meshed += 1;
                    self.report.triangles += idx.len() / 3;
                    self.report.vertices += pos.len();
                    self.out.push(FaceMesh {
                        face_id: *fid,
                        solid_id: f.solid,
                        positions: pos
                            .iter()
                            .flat_map(|p| [p.x as f32, p.y as f32, p.z as f32])
                            .collect(),
                        normals: nor
                            .iter()
                            .flat_map(|n| [n.x as f32, n.y as f32, n.z as f32])
                            .collect(),
                        indices: idx,
                    });
                }
                Err(why) => self.fail(why),
            }
        }

        if self.report.faces_meshed > 0 {
            self.report.bbox_min = lo;
            self.report.bbox_max = hi;
        }
        // Order matters: an opening is only recognisable once the faces round
        // it agree which way is out, and the patches themselves then have to be
        // turned to match, so the repair runs either side of the capping.
        self.report.rewound_triangles = orient_solid(&mut self.out);
        self.report.capped_holes = close_small_holes(&mut self.out, lo, hi);
        self.report.rewound_triangles += orient_solid(&mut self.out);
        let (open, bad) = mesh_defects(&self.out);
        self.report.open_edges = open;
        self.report.inconsistent_edges = bad;
    }

    /// Points along one `EDGE_CURVE`, from its start vertex to its end vertex,
    /// computed once and cached.
    fn edge_polyline(&mut self, edge: u32) -> Vec<Point3> {
        if let Some(p) = self.edge_pts.get(&edge) {
            return p.clone();
        }
        let (v1, v2, cid, same_sense) = match self.edges.get(&edge) {
            Some(e) => *e,
            None => return Vec::new(),
        };
        let pts = match self.curves.get(&cid) {
            Some(c) => sample_edge(c, v1, v2, same_sense, self.opts.chord_tol),
            None => {
                // No usable curve: a straight segment between the vertices is
                // the honest minimum, and the missing type is reported.
                vec![v1, v2]
            }
        };
        self.edge_pts.insert(edge, pts.clone());
        pts
    }

    fn finish(mut self) -> (Vec<FaceMesh>, ImportReport) {
        self.report.stray_loops = STRAY_LOOPS.with(|c| c.replace(0));
        if self.report.stray_loops > 0 {
            self.report.warnings.push(format!(
                "{} boundary loop(s) sat outside the face they were attached to and were ignored",
                self.report.stray_loops
            ));
        }
        let mut fails: Vec<(String, usize)> = self.fail_reasons.into_iter().collect();
        fails.sort_by(|a, b| b.1.cmp(&a.1));
        self.report.failures = fails
            .into_iter()
            .map(|(r, n)| if n > 1 { format!("{r} ({n} faces)") } else { r })
            .collect();

        let mut warns: Vec<(String, usize)> = self.warn_set.into_iter().collect();
        warns.sort_by(|a, b| b.1.cmp(&a.1));
        self.report.warnings.extend(
            warns
                .into_iter()
                .map(|(r, n)| if n > 1 { format!("{r} (x{n})") } else { r }),
        );

        (self.out, self.report)
    }
}

// ──────────────────────────── edge discretisation ──────────────────────────

/// Sample one edge between its two vertices.
///
/// The arc actually taken is decided here, and this is the step the old reader
/// skipped entirely: it stamped a full turn on every circular edge, so a 30°
/// fillet arc came back as a complete circle.
fn sample_edge(curve: &Curve, v1: Point3, v2: Point3, same_sense: bool, tol: f64) -> Vec<Point3> {
    let t1 = curve.param_of(v1);
    let t2 = curve.param_of(v2);

    let (a, b) = match curve.period() {
        Some(p) => {
            // Closed curve: the two vertices split it into two arcs, and the
            // edge's sense says which one this is.
            let closed = (v1 - v2).length() < 1e-9;
            if same_sense {
                let mut d = (t2 - t1).rem_euclid(p);
                if closed || d < 1e-12 {
                    d = p;
                }
                (t1, t1 + d)
            } else {
                let mut d = (t1 - t2).rem_euclid(p);
                if closed || d < 1e-12 {
                    d = p;
                }
                (t1, t1 - d)
            }
        }
        None => (t1, t2),
    };

    let mut pts = vec![curve.eval(a)];
    subdivide(curve, a, b, tol, 0, &mut pts);
    // Snap the ends to the vertices so neighbouring faces meet exactly, whatever
    // small error the parameter solve carries.
    if let Some(first) = pts.first_mut() {
        *first = v1;
    }
    if let Some(last) = pts.last_mut() {
        *last = v2;
    }
    pts
}

/// Recursive chord-tolerance sampling, appending everything after `a`.
fn subdivide(curve: &Curve, a: f64, b: f64, tol: f64, depth: u32, out: &mut Vec<Point3>) {
    let pa = curve.eval(a);
    let pb = curve.eval(b);
    let m = 0.5 * (a + b);
    let pm = curve.eval(m);
    let chord_mid = pa.lerp(pb, 0.5);
    let sag = (pm - chord_mid).length();
    // A closed span has a degenerate chord, so the sag test alone would accept
    // a full circle as a single segment; the angular cap stops that.
    let too_long = matches!(curve.period(), Some(p) if (b - a).abs() > p / 8.0);
    if depth < 9 && (sag > tol || too_long) {
        subdivide(curve, a, m, tol, depth + 1, out);
        subdivide(curve, m, b, tol, depth + 1, out);
    } else {
        out.push(pb);
    }
}

// ─────────────────────────────── face meshing ──────────────────────────────

type MeshOut = (Vec<Point3>, Vec<Vec3>, Vec<u32>);

/// One loop, carried in parameter space alongside its original 3-D points.
struct Ring {
    /// Unwrapped `(u, v)`.
    uv: Vec<[f64; 2]>,
    pts: Vec<Point3>,
    /// Turns of the periodic parameter over the whole loop: 0, +1 or -1.
    wrap_u: i32,
    wrap_v: i32,
}

fn mesh_face(surf: &Surf, rings3: &[Vec<Point3>], same_sense: bool, opts: &ImportOptions) -> Result<MeshOut, String> {
    let (pu, pv) = surf.periods();

    if rings3.is_empty() {
        // An untrimmed periodic surface — a whole torus or sphere — is the one
        // case where "draw the entire primitive" is the correct reading.
        return match (pu, pv) {
            (Some(_), Some(_)) => Ok(full_grid(surf, 0.0, TAU, -PI, PI, same_sense, opts)),
            _ => Err("face has no boundary and the surface is unbounded".into()),
        };
    }

    let mut rings: Vec<Ring> = Vec::with_capacity(rings3.len());
    for r in rings3 {
        rings.push(to_ring(surf, r, pu, pv));
    }

    let wraps_u = rings.iter().filter(|r| r.wrap_u != 0).count();
    let wraps_v = rings.iter().filter(|r| r.wrap_v != 0).count();

    // A band: two boundaries that each go once round the surface, with the
    // material between them. Cylinders and fillets are almost all of this.
    let mut axis = if wraps_u == 2 && wraps_v == 0 {
        Some(0usize)
    } else if wraps_v == 2 && wraps_u == 0 {
        Some(1usize)
    } else {
        None
    };
    // A torus whose tube is fatter than the ring it is swept round — the tight
    // corner of a serpentine — folds through its own axis, and there the second
    // parameter is undefined for a moment. A boundary that crosses that moment
    // picks up a spurious turn in the *other* direction, so a face that is a
    // plain band round one of them looks like it wraps both. Judge it by the
    // direction both boundaries agree on and let the band code have it.
    if axis.is_none() && rings.len() == 2 {
        if rings.iter().all(|r| r.wrap_v != 0) {
            axis = Some(1);
        } else if rings.iter().all(|r| r.wrap_u != 0) {
            axis = Some(0);
        }
    }
    if let Some(axis) = axis {
        // With no holes the band is a clean quad strip, which is both cheaper
        // and better shaped than anything a triangulator would produce.
        let period = if axis == 0 { pu } else { pv }.unwrap_or(TAU);
        if rings.len() == 2 {
            if let Ok(m) = band(surf, &rings, axis, period, same_sense, opts) {
                return Ok(m);
            }
        }
        // A filament that other filaments pass through has its side face
        // punched full of holes. Cut the band open at a seam and triangulate it
        // as an ordinary region with holes.
        match band_holes(surf, &rings, axis, period, same_sense, opts) {
            Ok(m) => return Ok(m),
            Err(e) => {
                if rings.len() == 2 {
                    // Last resort for a band whose boundaries wander: pair them
                    // up along their own length and loft. See `band_loft`.
                    return band_loft(surf, &rings, axis, period, same_sense, opts).map_err(|_| e);
                }
                return Err(e);
            }
        }
    }
    if wraps_u == 1 && rings.len() == 1 {
        if let Some(m) = cap_at_pole(surf, &rings[0], 0, pu.unwrap_or(TAU), same_sense, opts) {
            return m;
        }
    }
    if wraps_u > 0 || wraps_v > 0 {
        if std::env::var("STEP_IMPORT_DUMP").is_ok() {
            eprintln!(
                "[wrap] face #{} rings={} wraps_u={} wraps_v={}",
                DUMP_FACE.with(|c| c.get()),
                rings.len(),
                wraps_u,
                wraps_v
            );
            for (i, r) in rings.iter().enumerate() {
                let bb = |k: usize| {
                    let lo = r.uv.iter().map(|q| q[k]).fold(f64::MAX, f64::min);
                    let hi = r.uv.iter().map(|q| q[k]).fold(f64::MIN, f64::max);
                    (lo, hi)
                };
                let (u0, u1) = bb(0);
                let (v0, v1) = bb(1);
                let z0 = r.pts.iter().map(|p| p.z).fold(f64::MAX, f64::min);
                let z1 = r.pts.iter().map(|p| p.z).fold(f64::MIN, f64::max);
                eprintln!(
                    "        ring{i} n={} wrap=({},{}) u[{u0:.3},{u1:.3}] v[{v0:.3},{v1:.3}] z[{z0:.3},{z1:.3}]",
                    r.uv.len(),
                    r.wrap_u,
                    r.wrap_v
                );
            }
        }
        return Err(format!(
            "boundary wraps the surface {} time(s) in a pattern the importer cannot close",
            wraps_u + wraps_v
        ));
    }

    // Ordinary trimmed patch: the outer boundary is whichever loop encloses the
    // most area; the rest are holes.
    let (mu, mv) = surf.metric();
    let mut order: Vec<usize> = (0..rings.len()).collect();
    order.sort_by(|&a, &b| {
        let aa = ring_area(&rings[a].uv, mu, mv).abs();
        let bb = ring_area(&rings[b].uv, mu, mv).abs();
        bb.partial_cmp(&aa).unwrap_or(std::cmp::Ordering::Equal)
    });

    // Holes may have been unwrapped into a different period than the outer loop.
    let outer_idx = order[0];
    let base_u = mean(&rings[outer_idx].uv, 0);
    let base_v = mean(&rings[outer_idx].uv, 1);
    for (k, r) in rings.iter_mut().enumerate() {
        if k == outer_idx {
            continue;
        }
        if let Some(p) = pu {
            let d = p * ((base_u - mean(&r.uv, 0)) / p).round();
            for q in &mut r.uv {
                q[0] += d;
            }
        }
        if let Some(p) = pv {
            let d = p * ((base_v - mean(&r.uv, 1)) / p).round();
            for q in &mut r.uv {
                q[1] += d;
            }
        }
    }

    patch(surf, &rings, &order, same_sense, opts)
}

fn mean(uv: &[[f64; 2]], k: usize) -> f64 {
    if uv.is_empty() {
        0.0
    } else {
        uv.iter().map(|p| p[k]).sum::<f64>() / uv.len() as f64
    }
}

/// Project a 3-D loop into parameter space, unwrapping periodic jumps so the
/// loop reads as a continuous curve rather than one that teleports at the seam.
fn to_ring(surf: &Surf, pts: &[Point3], pu: Option<f64>, pv: Option<f64>) -> Ring {
    let n = pts.len();
    let mut uv = Vec::with_capacity(n);
    for p in pts {
        let (u, v) = surf.uv_of(*p);
        uv.push([u, v]);
    }
    let unwrap = |uv: &mut Vec<[f64; 2]>, k: usize, period: f64| -> i32 {
        for i in 1..uv.len() {
            let d = uv[i][k] - uv[i - 1][k];
            uv[i][k] -= period * (d / period).round();
        }
        // Closing step: how far the loop has travelled once it returns to its
        // own first point.
        let d = uv[0][k] - uv[uv.len() - 1][k];
        let closed = uv[uv.len() - 1][k] + d - period * (d / period).round();
        ((closed - uv[0][k]) / period).round() as i32
    };
    let wrap_u = pu.map(|p| unwrap(&mut uv, 0, p)).unwrap_or(0);
    let wrap_v = pv.map(|p| unwrap(&mut uv, 1, p)).unwrap_or(0);
    Ring {
        uv,
        pts: pts.to_vec(),
        wrap_u,
        wrap_v,
    }
}

fn ring_area(uv: &[[f64; 2]], mu: f64, mv: f64) -> f64 {
    let n = uv.len();
    let mut s = 0.0;
    for i in 0..n {
        let a = uv[i];
        let b = uv[(i + 1) % n];
        s += (a[0] * mu) * (b[1] * mv) - (b[0] * mu) * (a[1] * mv);
    }
    s * 0.5
}

// ── band: two boundaries that each wrap once ───────────────────────────────

/// Mesh the material between two boundaries that each go once round the
/// periodic direction `axis` (0 = u, 1 = v).
///
/// The two boundaries keep **exactly the points they arrived with**. An earlier
/// version resampled both onto a shared set of columns, which put a vertex in
/// the middle of an edge the neighbouring face had drawn as a single segment.
/// Nothing looks wrong — the point sits on the segment — but the two faces are
/// no longer joined along that edge, and the solid stops being closed. Instead
/// the two rims are stitched by walking them together and advancing whichever
/// one is behind, which needs no new points at all.
fn band(surf: &Surf, rings: &[Ring], axis: usize, period: f64, same_sense: bool, opts: &ImportOptions) -> Result<MeshOut, String> {
    let other = 1 - axis;
    let a = &rings[0];
    let b = &rings[1];
    let (lo_r, hi_r) = if mean(&a.uv, other) <= mean(&b.uv, other) {
        (a, b)
    } else {
        (b, a)
    };

    let lo = monotone_columns(lo_r, axis, period).ok_or("boundary doubles back over the seam")?;
    let mut hi = monotone_columns(hi_r, axis, period).ok_or("boundary doubles back over the seam")?;
    let d = window_shift(hi[0].0, lo[0].0, period);
    for c in &mut hi {
        c.0 += d;
    }
    if lo.len() < 3 || hi.len() < 3 {
        return Err("band boundary has too few points".into());
    }

    let probe = lo[lo.len() / 2].0;
    let rows = rows_needed(surf, axis, probe, interp_at(&lo, probe, period), interp_at(&hi, probe, period), opts).max(1);
    let merged_len = lo.len() + hi.len();
    let rows = if merged_len * rows > opts.max_triangles_per_face {
        (opts.max_triangles_per_face / merged_len.max(1)).max(1)
    } else {
        rows
    };

    let uv_of = |c: f64, o: f64| if axis == 0 { (c, o) } else { (o, c) };
    let mut pos: Vec<Point3> = Vec::new();
    let mut nor: Vec<Vec3> = Vec::new();
    let push = |p: Point3, c: f64, o: f64, pos: &mut Vec<Point3>, nor: &mut Vec<Vec3>| -> u32 {
        let (u, v) = uv_of(c, o);
        pos.push(p);
        nor.push(surf.normal(u, v));
        (pos.len() - 1) as u32
    };

    // The two rims, as they came in.
    let row_lo: Vec<(f64, u32)> = lo
        .iter()
        .map(|c| (c.0, push(c.2, c.0, c.1, &mut pos, &mut nor)))
        .collect();
    let row_hi: Vec<(f64, u32)> = hi
        .iter()
        .map(|c| (c.0, push(c.2, c.0, c.1, &mut pos, &mut nor)))
        .collect();

    // Interior rows, where the surface bends across the band. These are the
    // face's own points, so they may sit wherever the tolerance wants them.
    let mut inner: Vec<Vec<(f64, u32)>> = Vec::new();
    if rows > 1 {
        // Both rims are folded into one turn before merging. The upper rim was
        // shifted to start inside the lower rim's turn, so its own columns can
        // run past the end of it; left unfolded the interior rows would wind
        // twice round the surface and the strip would tie itself in knots.
        let base = lo[0].0;
        let mut cols: Vec<f64> = lo
            .iter()
            .chain(hi.iter())
            .map(|c| base + (c.0 - base).rem_euclid(period))
            .collect();
        cols.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        cols.dedup_by(|x, y| (*x - *y).abs() < 1e-9);
        for k in 1..rows {
            let t = k as f64 / rows as f64;
            let mut row = Vec::with_capacity(cols.len());
            for &c in &cols {
                let o0 = interp_at(&lo, c, period);
                let o1 = interp_at(&hi, c, period);
                let o = o0 + (o1 - o0) * t;
                let (u, v) = uv_of(c, o);
                row.push((c, push(surf.eval(u, v), c, o, &mut pos, &mut nor)));
            }
            inner.push(row);
        }
    }

    let mut idx: Vec<u32> = Vec::new();
    // A triangle wound counter-clockwise in (column, other) is counter-clockwise
    // in (u, v) only when the columns run along u.
    let keep = (axis == 0) == same_sense;
    let mut tri = |a: u32, b: u32, c: u32, idx: &mut Vec<u32>| {
        if keep {
            idx.extend([a, b, c]);
        } else {
            idx.extend([a, c, b]);
        }
    };

    let mut prev = &row_lo;
    for row in inner.iter() {
        stitch(prev, row, period, &mut idx, &mut tri);
        prev = row;
    }
    stitch(prev, &row_hi, period, &mut idx, &mut tri);

    if !same_sense {
        for n in &mut nor {
            *n = -*n;
        }
    }
    Ok((pos, nor, idx))
}

/// Value of a boundary at column `c`, interpolated along the polyline.
fn interp_at(cols: &[(f64, f64, Point3)], c: f64, period: f64) -> f64 {
    let base = cols[0].0;
    let x = base + (c - base).rem_euclid(period);
    let mut i = 0;
    while i + 2 < cols.len() && cols[i + 1].0 <= x {
        i += 1;
    }
    let (x0, q0, _) = cols[i];
    let (x1, q1, _) = cols[i + 1];
    let d = x1 - x0;
    let t = if d.abs() < 1e-15 { 0.0 } else { ((x - x0) / d).clamp(0.0, 1.0) };
    q0 + (q1 - q0) * t
}

/// Triangles between two closed runs that each go once round the band.
///
/// Walks both at once and always advances the one lying further back, so every
/// point of both runs is used and none is invented.
///
/// The two runs are first brought onto the same turn and rotated to start at
/// the same place. Skipping that looks harmless — the walk still consumes both
/// runs and still produces a closed strip — but it pairs each point of one rim
/// with a point of the other some way round the tube, and the side comes out
/// **twisted**: a fan at each end, a pinched waist, and a solid that has lost a
/// fifth of its volume while looking very nearly right.
fn stitch(
    a: &[(f64, u32)],
    b: &[(f64, u32)],
    period: f64,
    idx: &mut Vec<u32>,
    tri: &mut impl FnMut(u32, u32, u32, &mut Vec<u32>),
) {
    let (na, nb) = (a.len(), b.len());
    if na < 2 || nb < 2 {
        return;
    }
    let base = a[0].0;
    // `b` folded into `a`'s turn and rotated to begin where `a` begins.
    let mut bb: Vec<(f64, u32)> = b
        .iter()
        .map(|&(c, i)| (base + (c - base).rem_euclid(period), i))
        .collect();
    bb.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));

    let col_a = |i: usize| a[i % na].0 + if i >= na { period } else { 0.0 };
    let col_b = |j: usize| bb[j % nb].0 + if j >= nb { period } else { 0.0 };

    let (mut i, mut j) = (0usize, 0usize);
    while i < na || j < nb {
        let advance_a = if i >= na {
            false
        } else if j >= nb {
            true
        } else {
            col_a(i + 1) <= col_b(j + 1)
        };
        if advance_a {
            tri(a[i % na].1, a[(i + 1) % na].1, bb[j % nb].1, idx);
            i += 1;
        } else {
            tri(a[i % na].1, bb[(j + 1) % nb].1, bb[j % nb].1, idx);
            j += 1;
        }
    }
}

/// A boundary expressed as `(column, other, point)`, ordered so the column
/// advances once round the surface. `None` when the boundary is not
/// single-valued in the column, which the quad-strip construction cannot
/// represent — the seam-cut path takes over in that case.
fn monotone_columns(r: &Ring, axis: usize, period: f64) -> Option<Vec<(f64, f64, Point3)>> {
    let n = r.uv.len();
    if n < 3 {
        return None;
    }
    let mut v: Vec<(f64, f64, Point3)> = (0..n)
        .map(|i| (r.uv[i][axis], r.uv[i][1 - axis], r.pts[i]))
        .collect();
    if v[n - 1].0 < v[0].0 {
        v.reverse();
    }
    if !(v[n - 1].0 - v[0].0 > period * 0.5) {
        return None;
    }
    for i in 1..n {
        if v[i].0 < v[i - 1].0 - 1e-12 {
            return None;
        }
    }
    Some(v)
}

/// Emit one grid cell with the winding that makes its 3-D normal point out of
/// the solid.
///
/// A triangle wound counter-clockwise in `(u, v)` carries the surface's own
/// normal, `du x dv`; STEP's `same_sense` flag says whether that is the outward
/// direction or its opposite. Which index order is counter-clockwise depends on
/// whether the grid's columns run along `u` or along `v` — getting that
/// backwards turns a solid inside out without moving a single vertex, which is
/// invisible until something measures the volume.
fn push_quad(idx: &mut Vec<u32>, a0: u32, a1: u32, b0: u32, b1: u32, axis: usize, same_sense: bool) {
    // (a0, a1) steps along the row, (a0, b0) along the column.
    let (t0, t1) = if axis == 0 {
        ([a0, b0, b1], [a0, b1, a1])
    } else {
        ([a0, a1, b1], [a0, b1, b0])
    };
    if same_sense {
        idx.extend(t0);
        idx.extend(t1);
    } else {
        idx.extend([t0[0], t0[2], t0[1]]);
        idx.extend([t1[0], t1[2], t1[1]]);
    }
}

/// Loft between two boundaries that each go once round, pairing them by their
/// own length rather than by the surface parameter.
///
/// This is the fallback for a band whose boundaries double back on themselves,
/// where there is no single column at which to lay a straight cut. It is an
/// **approximation**: the boundaries are reproduced exactly, and every interior
/// point is evaluated on the surface, but the way the two are paired up is a
/// choice this code makes rather than something the file states. It runs only
/// after both exact routes have failed, and the faces it draws are counted, so
/// a model that leans on it says so.
fn band_loft(surf: &Surf, rings: &[Ring], axis: usize, period: f64, same_sense: bool, opts: &ImportOptions) -> Result<MeshOut, String> {
    let (mu, mv) = surf.metric();
    let scale = |p: [f64; 2]| [p[0] * mu, p[1] * mv];

    // Cumulative length round each boundary, normalised to one turn.
    let arc = |r: &Ring| -> Vec<f64> {
        let n = r.uv.len();
        let mut t = Vec::with_capacity(n + 1);
        t.push(0.0);
        for i in 0..n {
            let a = scale(r.uv[i]);
            let b = scale(r.uv[(i + 1) % n]);
            let d = (a[0] - b[0]).hypot(a[1] - b[1]);
            let last = *t.last().unwrap();
            t.push(last + d);
        }
        let total = *t.last().unwrap();
        if total <= 0.0 {
            return vec![0.0; n + 1];
        }
        t.iter().map(|x| x / total).collect()
    };

    let a = &rings[0];
    let b = &rings[1];
    let (ta, tb) = (arc(a), arc(b));
    if a.uv.len() < 3 || b.uv.len() < 3 {
        return Err("band boundary has too few points".into());
    }

    // Start the second boundary at whichever of its points is nearest the first
    // one's start, so the loft does not begin with a twist.
    let start_b = (0..b.pts.len())
        .min_by(|&i, &j| {
            (b.pts[i] - a.pts[0])
                .length()
                .partial_cmp(&(b.pts[j] - a.pts[0]).length())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap_or(0);
    let off = tb[start_b];

    // Bring the second boundary onto the same turn of the surface.
    let shift = window_shift(mean(&b.uv, axis), mean(&a.uv, axis) - period * 0.5, period);

    let uv_a = |k: usize| a.uv[k % a.uv.len()];
    let uv_b = |k: usize| {
        let mut p = b.uv[k % b.uv.len()];
        p[axis] += shift;
        p
    };

    let rows = {
        let ka = 0usize;
        let kb = start_b;
        let (pa, pb) = (uv_a(ka), uv_b(kb));
        let mut n = 1usize;
        for cand in [1usize, 2, 3, 4, 6, 8, 12, 16, 24, 32] {
            let mut worst: f64 = 0.0;
            for k in 0..cand {
                let f0 = k as f64 / cand as f64;
                let f1 = (k + 1) as f64 / cand as f64;
                let at = |f: f64| surf.eval(pa[0] + (pb[0] - pa[0]) * f, pa[1] + (pb[1] - pa[1]) * f);
                let m = at(0.5 * (f0 + f1));
                worst = worst.max((m - at(f0).lerp(at(f1), 0.5)).length());
            }
            n = cand;
            if worst <= opts.chord_tol {
                break;
            }
        }
        n
    };

    let mut pos: Vec<Point3> = Vec::new();
    let mut nor: Vec<Vec3> = Vec::new();
    let push = |p: Point3, uv: [f64; 2], pos: &mut Vec<Point3>, nor: &mut Vec<Vec3>| -> u32 {
        pos.push(p);
        nor.push(surf.normal(uv[0], uv[1]));
        (pos.len() - 1) as u32
    };

    let row_a: Vec<(f64, u32)> = (0..a.uv.len())
        .map(|i| (ta[i], push(a.pts[i], uv_a(i), &mut pos, &mut nor)))
        .collect();
    let row_b: Vec<(f64, u32)> = (0..b.uv.len())
        .map(|i| {
            let k = (start_b + i) % b.uv.len();
            let t = (tb[k] - off).rem_euclid(1.0);
            (t, push(b.pts[k], uv_b(k), &mut pos, &mut nor))
        })
        .collect();

    // Interior rows follow the surface between the paired boundary points.
    let mut inner: Vec<Vec<(f64, u32)>> = Vec::new();
    if rows > 1 {
        let mut ts: Vec<f64> = row_a.iter().chain(row_b.iter()).map(|c| c.0).collect();
        ts.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        ts.dedup_by(|x, y| (*x - *y).abs() < 1e-12);
        let at = |seq: &[f64], uv: &dyn Fn(usize) -> [f64; 2], t: f64| -> [f64; 2] {
            let n = seq.len();
            let mut i = 0;
            while i + 1 < n && seq[i + 1] <= t {
                i += 1;
            }
            let (t0, t1) = (seq[i], if i + 1 < n { seq[i + 1] } else { seq[i] + 1.0 });
            let f = if (t1 - t0).abs() < 1e-15 { 0.0 } else { (t - t0) / (t1 - t0) };
            let (p0, p1) = (uv(i), uv(i + 1));
            [p0[0] + (p1[0] - p0[0]) * f, p0[1] + (p1[1] - p0[1]) * f]
        };
        let ts_a: Vec<f64> = row_a.iter().map(|c| c.0).collect();
        let ts_b: Vec<f64> = row_b.iter().map(|c| c.0).collect();
        let uvb = |i: usize| uv_b(start_b + i);
        for k in 1..rows {
            let f = k as f64 / rows as f64;
            let mut row = Vec::with_capacity(ts.len());
            for &t in &ts {
                let pa = at(&ts_a, &uv_a, t);
                let pb = at(&ts_b, &uvb, t);
                let uv = [pa[0] + (pb[0] - pa[0]) * f, pa[1] + (pb[1] - pa[1]) * f];
                row.push((t, push(surf.eval(uv[0], uv[1]), uv, &mut pos, &mut nor)));
            }
            inner.push(row);
        }
    }

    let keep = (axis == 0) == same_sense;
    let mut tri = |x: u32, y: u32, z: u32, idx: &mut Vec<u32>| {
        if keep {
            idx.extend([x, y, z]);
        } else {
            idx.extend([x, z, y]);
        }
    };
    let mut idx: Vec<u32> = Vec::new();
    let mut prev = &row_a;
    for row in inner.iter() {
        stitch(prev, row, 1.0, &mut idx, &mut tri);
        prev = row;
    }
    stitch(prev, &row_b, 1.0, &mut idx, &mut tri);

    if !same_sense {
        for n in &mut nor {
            *n = -*n;
        }
    }
    Ok((pos, nor, idx))
}

/// A loop as `(column, other, point)`, ordered so the column advances.
fn ring_seq(r: &Ring, axis: usize) -> Vec<(f64, f64, Point3)> {
    let n = r.uv.len();
    let mut v: Vec<(f64, f64, Point3)> = (0..n)
        .map(|i| (r.uv[i][axis], r.uv[i][1 - axis], r.pts[i]))
        .collect();
    if n >= 2 && v[n - 1].0 < v[0].0 {
        v.reverse();
    }
    v
}

/// How many steps across the band keep the surface within tolerance.
fn rows_needed(surf: &Surf, axis: usize, c: f64, q0: f64, q1: f64, opts: &ImportOptions) -> usize {
    if !surf.is_curved() {
        return 1;
    }
    let at = |q: f64| {
        if axis == 0 {
            surf.eval(c, q)
        } else {
            surf.eval(q, c)
        }
    };
    for n in [1usize, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64] {
        let mut worst: f64 = 0.0;
        for k in 0..n {
            let a = q0 + (q1 - q0) * k as f64 / n as f64;
            let b = q0 + (q1 - q0) * (k + 1) as f64 / n as f64;
            let m = at(0.5 * (a + b));
            worst = worst.max((m - at(a).lerp(at(b), 0.5)).length());
        }
        if worst <= opts.chord_tol {
            return n;
        }
    }
    64
}

/// A band with holes punched through it — a filament that other filaments pass
/// through, which is most of a fused scaffold.
///
/// A band is a tube, not a disk, so no polygon triangulator can take it
/// directly. It is cut open along a seam that runs from one boundary to the
/// other, which turns the tube into a disk with holes. The cut is placed at
/// **existing sample points on both boundaries**, never at an interpolated
/// position, so no vertex is invented on an edge that a neighbouring face has
/// already drawn — that is what keeps the two faces closed against each other.
/// Its column is chosen as far as possible from every hole, and if the outline
/// that results crosses itself, the next-best column is tried instead.
fn band_holes(surf: &Surf, rings: &[Ring], axis: usize, period: f64, same_sense: bool, opts: &ImportOptions) -> Result<MeshOut, String> {
    let wrap_of = |r: &Ring| if axis == 0 { r.wrap_u } else { r.wrap_v };
    let rims: Vec<usize> = (0..rings.len()).filter(|&i| wrap_of(&rings[i]) != 0).collect();
    if rims.len() != 2 {
        return Err(format!(
            "band has {} boundaries going round it, and only two can be closed",
            rims.len()
        ));
    }

    let s0 = ring_seq(&rings[rims[0]], axis);
    let mut s1 = ring_seq(&rings[rims[1]], axis);
    let base = s0[0].0;
    // Both rims must be read on the *same* turn of the surface. Lining them up
    // by their average column is not enough: rounding can land a whole period
    // out, and the band then unrolls as two turns laid end to end instead of
    // one, which no triangulation can make sense of.
    let d1 = window_shift(s1[0].0, base, period);
    for c in &mut s1 {
        c.0 += d1;
    }

    // Where the holes sit, so the cut can avoid them.
    let mut holes: Vec<Vec<(f64, f64, Point3)>> = Vec::new();
    for i in 0..rings.len() {
        if rims.contains(&i) {
            continue;
        }
        let mut h = ring_seq(&rings[i], axis);
        let d = window_shift(mean_col(&h), base, period);
        for c in &mut h {
            c.0 += d;
        }
        holes.push(h);
    }

    let span = |h: &[(f64, f64, Point3)]| {
        let lo = h.iter().map(|c| c.0).fold(f64::INFINITY, f64::min);
        let hi = h.iter().map(|c| c.0).fold(f64::NEG_INFINITY, f64::max);
        (lo, hi)
    };
    let spans: Vec<(f64, f64)> = holes.iter().map(|h| span(h)).collect();
    let circ = |x: f64| {
        let t = x.rem_euclid(period);
        t.min(period - t)
    };
    let mut ranked: Vec<(f64, f64)> = (0..256)
        .map(|i| {
            let c = base + period * i as f64 / 256.0;
            let mut d = f64::INFINITY;
            for &(lo, hi) in &spans {
                let inside = (c - lo).rem_euclid(period) <= (hi - lo).rem_euclid(period);
                d = d.min(if inside { -1.0 } else { circ(c - lo).min(circ(c - hi)) });
            }
            (d, c)
        })
        .collect();
    ranked.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

    let near = |s: &[(f64, f64, Point3)], c: f64| -> usize {
        (0..s.len())
            .min_by(|&a, &b| {
                circ(s[a].0 - c)
                    .partial_cmp(&circ(s[b].0 - c))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(0)
    };
    let to_uv = |c: &(f64, f64, Point3)| -> ([f64; 2], Point3) {
        if axis == 0 {
            ([c.0, c.1], c.2)
        } else {
            ([c.1, c.0], c.2)
        }
    };

    let mut chosen: Option<Vec<([f64; 2], Point3)>> = None;
    let mut holes_used: Vec<Vec<(f64, f64, Point3)>> = Vec::new();
    for &(_, cut) in ranked.iter().take(24) {
        let mut a_fwd = rotate_forward(&s0, near(&s0, cut), period);
        let da = window_shift(a_fwd[0].0, cut - period * 0.5, period);
        for c in &mut a_fwd {
            c.0 += da;
        }
        let c0 = a_fwd[0].0;

        let mut b_rev = rotate_forward(&s1, near(&s1, cut), period);
        let db = window_shift(b_rev[0].0, c0 - period * 0.5, period);
        for c in &mut b_rev {
            c.0 += db;
        }
        b_rev.reverse();

        // The cut itself is subdivided rather than left as one long segment. Two
        // reasons: a single edge spanning the whole face gives the hole-bridging
        // nothing sensible to attach to, and the cut follows the surface only if
        // it is sampled. Both copies get the same fractions, so they stay
        // congruent a period apart and the face still closes against itself.
        let seam_n = 24usize;
        let (ca, oa) = (a_fwd[0].0, a_fwd[0].1);
        let (cb, ob) = (b_rev[b_rev.len() - 1].0, b_rev[b_rev.len() - 1].1);
        let seam_at = |k: usize, shift: f64| -> (f64, f64, Point3) {
            let t = k as f64 / seam_n as f64;
            let c = ca + (cb - ca) * t + shift;
            let o = oa + (ob - oa) * t;
            let (u, v) = if axis == 0 { (c, o) } else { (o, c) };
            (c, o, surf.eval(u, v))
        };

        let mut outer: Vec<([f64; 2], Point3)> = Vec::with_capacity(a_fwd.len() + b_rev.len() + 2 * seam_n);
        outer.extend(a_fwd.iter().map(to_uv));
        for k in 1..seam_n {
            outer.push(to_uv(&seam_at(k, period)));
        }
        outer.extend(b_rev.iter().map(to_uv));
        for k in (1..seam_n).rev() {
            outer.push(to_uv(&seam_at(k, 0.0)));
        }

        let shifted: Vec<Vec<(f64, f64, Point3)>> = holes
            .iter()
            .map(|h| {
                let d = window_shift(mean_col(h), c0, period);
                h.iter().map(|c| (c.0 + d, c.1, c.2)).collect()
            })
            .collect();

        let hole_uv: Vec<Vec<[f64; 2]>> = shifted
            .iter()
            .map(|h| h.iter().map(|c| to_uv(c).0).collect())
            .collect();
        let outer_uv: Vec<[f64; 2]> = outer.iter().map(|(p, _)| *p).collect();
        let ok = simple_outline(&outer_uv, &hole_uv);
        if ok || chosen.is_none() {
            chosen = Some(outer);
            holes_used = shifted;
        }
        if ok {
            break;
        }
    }
    let outer = chosen.ok_or("the band could not be cut open")?;

    let mut loops: Vec<Vec<([f64; 2], Point3)>> = vec![outer];
    for h in holes_used {
        loops.push(h.iter().map(to_uv).collect());
    }

    mesh_contours(surf, &loops, same_sense, opts)
}

/// Whole periods to add to `x` to bring it into `[lo, lo + period)`.
fn window_shift(x: f64, lo: f64, period: f64) -> f64 {
    lo + (x - lo).rem_euclid(period) - x
}

fn mean_col(s: &[(f64, f64, Point3)]) -> f64 {
    if s.is_empty() {
        0.0
    } else {
        s.iter().map(|c| c.0).sum::<f64>() / s.len() as f64
    }
}

/// One full turn of a closed boundary starting at sample `s`, with the columns
/// past the wrap-around advanced by a period so the run reads as monotone.
fn rotate_forward(seq: &[(f64, f64, Point3)], s: usize, period: f64) -> Vec<(f64, f64, Point3)> {
    let n = seq.len();
    let mut out = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let idx = (s + i) % n;
        let bump = if s + i >= n { period } else { 0.0 };
        out.push((seq[idx].0 + bump, seq[idx].1, seq[idx].2));
    }
    out
}

/// A boundary that wraps once with nothing on the far side: the surface closes
/// itself at a pole, as a sphere does at its axis or a cone at its apex.
fn cap_at_pole(surf: &Surf, r: &Ring, axis: usize, period: f64, same_sense: bool, opts: &ImportOptions) -> Option<Result<MeshOut, String>> {
    let poles = surf.poles();
    if poles.is_empty() {
        return None;
    }
    let q_mean = mean(&r.uv, 1 - axis);
    let pole = *poles
        .iter()
        .min_by(|a, b| {
            (**a - q_mean)
                .abs()
                .partial_cmp(&(**b - q_mean).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap();

    let cols = match monotone_columns(r, axis, period) {
        Some(c) => c,
        None => return Some(Err("boundary around the pole is not single valued".into())),
    };
    let rows = rows_needed(surf, axis, cols[cols.len() / 2].0, cols[cols.len() / 2].1, pole, opts).max(1);

    let mut pos = Vec::new();
    let mut nor = Vec::new();
    let nc = cols.len() - 1; // last column repeats the first
    let stride = rows + 1;
    for c in cols.iter().take(nc) {
        for k in 0..=rows {
            let t = k as f64 / rows as f64;
            let q = c.1 + (pole - c.1) * t;
            let uv = if axis == 0 { (c.0, q) } else { (q, c.0) };
            pos.push(if k == 0 { c.2 } else { surf.eval(uv.0, uv.1) });
            nor.push(surf.normal(uv.0, uv.1));
        }
    }
    let mut idx = Vec::new();
    for i in 0..nc {
        let i1 = (i + 1) % nc;
        for k in 0..rows {
            let a0 = (i * stride + k) as u32;
            let a1 = (i * stride + k + 1) as u32;
            let b0 = (i1 * stride + k) as u32;
            let b1 = (i1 * stride + k + 1) as u32;
            push_quad(&mut idx, a0, a1, b0, b1, axis, same_sense);
        }
    }
    if !same_sense {
        for n in &mut nor {
            *n = -*n;
        }
    }
    Some(Ok((pos, nor, idx)))
}

/// Whole periodic surface, used only when a face carries no boundary at all.
fn full_grid(surf: &Surf, u0: f64, u1: f64, v0: f64, v1: f64, same_sense: bool, opts: &ImportOptions) -> MeshOut {
    let (mu, mv) = surf.metric();
    let nu = (((u1 - u0) * mu / (8.0 * opts.chord_tol).max(1e-6)).ceil() as usize).clamp(8, 256);
    let nv = (((v1 - v0) * mv / (8.0 * opts.chord_tol).max(1e-6)).ceil() as usize).clamp(8, 256);
    let mut pos = Vec::new();
    let mut nor = Vec::new();
    for i in 0..nu {
        for j in 0..=nv {
            let u = u0 + (u1 - u0) * i as f64 / nu as f64;
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            pos.push(surf.eval(u, v));
            nor.push(if same_sense { surf.normal(u, v) } else { -surf.normal(u, v) });
        }
    }
    let mut idx = Vec::new();
    let stride = nv + 1;
    for i in 0..nu {
        let i1 = (i + 1) % nu;
        for j in 0..nv {
            let a0 = (i * stride + j) as u32;
            let a1 = (i * stride + j + 1) as u32;
            let b0 = (i1 * stride + j) as u32;
            let b1 = (i1 * stride + j + 1) as u32;
            push_quad(&mut idx, a0, a1, b0, b1, 0, same_sense);
        }
    }
    (pos, nor, idx)
}

// ── patch: an ordinary trimmed region ──────────────────────────────────────

fn patch(surf: &Surf, rings: &[Ring], order: &[usize], same_sense: bool, opts: &ImportOptions) -> Result<MeshOut, String> {
    let loops: Vec<Vec<([f64; 2], Point3)>> = order
        .iter()
        .map(|&k| {
            let r = &rings[k];
            (0..r.uv.len()).map(|i| (r.uv[i], r.pts[i])).collect()
        })
        .collect();
    mesh_contours(surf, &loops, same_sense, opts)
}

/// Triangulate a face given its boundary in parameter space, outer loop first.
///
/// Everything that is not a clean quad strip ends up here: planar caps with
/// holes, trimmed patches, and bands that have been cut open at a seam.
fn mesh_contours(
    surf: &Surf,
    loops: &[Vec<([f64; 2], Point3)>],
    same_sense: bool,
    opts: &ImportOptions,
) -> Result<MeshOut, String> {
    let (mu, mv) = surf.metric();

    // One shared point list; the triangulator works on indices into it.
    let mut uv: Vec<[f64; 2]> = Vec::new();
    let mut pos: Vec<Point3> = Vec::new();
    let mut contours: Vec<Vec<usize>> = Vec::new();
    for lp in loops {
        let mut c = Vec::with_capacity(lp.len());
        for (p, w) in lp {
            c.push(uv.len());
            uv.push([p[0] * mu, p[1] * mv]);
            pos.push(*w);
        }
        contours.push(c);
    }
    if contours.is_empty() || contours[0].len() < 3 {
        return Err("face outline has fewer than three points".into());
    }

    let outer = contours[0].clone();
    let mut holes: Vec<Vec<usize>> = contours[1..].to_vec();

    // Loops the file attached to this face that are not on it.
    //
    // A scaffold exporter hands the *straight cylinder* of a filament the
    // crossing window that belongs to the elbow beside it: the loop sits past
    // the face's own end circle, at a negative axial parameter. A hole outside
    // the outline has no reading — it encloses none of this face — and keeping
    // it leaves an outline no triangulator can close, which is why the entire
    // face used to be dropped. That is what emptied the section outline at
    // every crossing. Dropping the stray loop keeps the face.
    let strays = holes.len();
    holes.retain(|h| h.iter().any(|&i| point_in_polygon(&uv, &outer, uv[i])));
    let strays = strays - holes.len();
    if strays > 0 {
        STRAY_LOOPS.with(|c| c.set(c.get() + strays));
    }

    let mut tris = match triangulate(&uv, &outer, &holes)
        .or_else(|e| {
            // A window that straddles the face's own rim is not a hole at all —
            // it is a notch, and the file states it as a loop crossing the
            // outline. Triangulating the face without it leaves a patch of
            // surface where a neighbour passes through; that is a fraction of a
            // millimetre of skin inside another solid, which the section's
            // union absorbs and the render never shows. Losing the whole face
            // instead leaves a hole a millimetre wide in the outline.
            let kept: Vec<Vec<usize>> = holes
                .iter()
                .filter(|h| !loop_crosses(&uv, &outer, h))
                .cloned()
                .collect();
            if kept.len() == holes.len() {
                return Err(e);
            }
            STRAY_LOOPS.with(|c| c.set(c.get() + holes.len() - kept.len()));
            triangulate(&uv, &outer, &kept)
        }) {
        Ok(t) => t,
        Err(e) => {
            // A face the triangulator gave up on is worth looking at rather
            // than guessing about, so it can be dumped on request.
            if let Ok(dir) = std::env::var("STEP_IMPORT_DUMP") {
                let mut txt = String::new();
                for (li, c) in contours.iter().enumerate() {
                    for &i in c {
                        txt.push_str(&format!(
                            "{li} {:.9} {:.9} {:.6} {:.6} {:.6}
",
                            uv[i][0], uv[i][1], pos[i].x, pos[i].y, pos[i].z
                        ));
                    }
                }
                let n = DUMP_FACE.with(|c| c.get());
                let _unused = std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .map(|d| d.subsec_nanos())
                    .unwrap_or(0);
                let _ = std::fs::write(format!("{dir}/face_{n}.txt"), txt);
            }
            return Err(e);
        }
    };

    // Boundary edges must not be moved by refinement, or the faces next door
    // pull apart from this one.
    let mut boundary: HashSet<(usize, usize)> = HashSet::new();
    for c in &contours {
        for i in 0..c.len() {
            let a = c[i];
            let b = c[(i + 1) % c.len()];
            boundary.insert((a.min(b), a.max(b)));
        }
    }

    if surf.is_curved() {
        refine(surf, &mut uv, &mut pos, &mut tris, &mut boundary, mu, mv, opts);
    }

    let nor: Vec<Vec3> = uv
        .iter()
        .map(|p| {
            let n = surf.normal(p[0] / mu, p[1] / mv);
            if same_sense {
                n
            } else {
                -n
            }
        })
        .collect();

    let mut idx = Vec::with_capacity(tris.len() * 3);
    for t in tris {
        if same_sense {
            idx.extend([t[0] as u32, t[1] as u32, t[2] as u32]);
        } else {
            idx.extend([t[0] as u32, t[2] as u32, t[1] as u32]);
        }
    }
    Ok((pos, nor, idx))
}

/// Subdivide until the mesh follows the surface, **without ever moving the
/// boundary**.
///
/// A face may only add points in its own interior. The moment refinement splits
/// an edge the face shares with its neighbour, the neighbour — which refines on
/// its own schedule, or is a plane and does not refine at all — no longer has a
/// matching vertex there, and the two meshes stop being joined. The surfaces
/// still touch, so nothing looks wrong, but the solid is no longer closed and
/// anything that walks the mesh afterwards falls through the gap.
///
/// So the split is decided per *edge*, globally: an interior edge either splits
/// for both triangles that own it or for neither, and the four standard
/// templates cover 3, 2, 1 and 0 split edges. A triangle bounded entirely by
/// shared edges cannot split any of them, and is refined at its centroid, which
/// needs no agreement with anybody.
#[allow(clippy::too_many_arguments)]
fn refine(
    surf: &Surf,
    uv: &mut Vec<[f64; 2]>,
    pos: &mut Vec<Point3>,
    tris: &mut Vec<[usize; 3]>,
    boundary: &mut HashSet<(usize, usize)>,
    mu: f64,
    mv: f64,
    opts: &ImportOptions,
) {
    let key = |a: usize, b: usize| (a.min(b), a.max(b));

    for _ in 0..6 {
        // How far a triangle strays from the surface it stands for, measured at
        // its centroid and at the middle of each side.
        let sag = |uv: &Vec<[f64; 2]>, pos: &Vec<Point3>, t: &[usize; 3]| -> f64 {
            let mut worst: f64 = 0.0;
            let probe = |a: usize, b: usize, c: usize, w: f64| {
                let p = [
                    (uv[a][0] + uv[b][0] + uv[c][0]) * w,
                    (uv[a][1] + uv[b][1] + uv[c][1]) * w,
                ];
                let q = Point3::new(
                    (pos[a].x + pos[b].x + pos[c].x) * w,
                    (pos[a].y + pos[b].y + pos[c].y) * w,
                    (pos[a].z + pos[b].z + pos[c].z) * w,
                );
                (surf.eval(p[0] / mu, p[1] / mv) - q).length()
            };
            worst = worst.max(probe(t[0], t[1], t[2], 1.0 / 3.0));
            for e in 0..3 {
                let (a, b) = (t[e], t[(e + 1) % 3]);
                worst = worst.max(probe(a, b, b, 0.5));
            }
            worst
        };

        let need: Vec<bool> = tris.iter().map(|t| sag(uv, pos, t) > opts.chord_tol).collect();
        if !need.iter().any(|&x| x) || tris.len() * 4 > opts.max_triangles_per_face {
            return;
        }

        // Marking is per edge, so neighbours cannot disagree.
        let mut marked: HashSet<(usize, usize)> = HashSet::new();
        for (i, t) in tris.iter().enumerate() {
            if !need[i] {
                continue;
            }
            for e in 0..3 {
                let k = key(t[e], t[(e + 1) % 3]);
                if !boundary.contains(&k) {
                    marked.insert(k);
                }
            }
        }

        let mut mid: HashMap<(usize, usize), usize> = HashMap::new();
        let mut split = |uv: &mut Vec<[f64; 2]>, pos: &mut Vec<Point3>, a: usize, b: usize| -> usize {
            *mid.entry(key(a, b)).or_insert_with(|| {
                let p = [0.5 * (uv[a][0] + uv[b][0]), 0.5 * (uv[a][1] + uv[b][1])];
                let w = surf.eval(p[0] / mu, p[1] / mv);
                uv.push(p);
                pos.push(w);
                uv.len() - 1
            })
        };

        let mut out: Vec<[usize; 3]> = Vec::with_capacity(tris.len() * 2);
        for (i, t) in tris.iter().enumerate() {
            let m: [bool; 3] = [
                marked.contains(&key(t[0], t[1])),
                marked.contains(&key(t[1], t[2])),
                marked.contains(&key(t[2], t[0])),
            ];
            let count = m.iter().filter(|&&x| x).count();
            match count {
                3 => {
                    let m0 = split(uv, pos, t[0], t[1]);
                    let m1 = split(uv, pos, t[1], t[2]);
                    let m2 = split(uv, pos, t[2], t[0]);
                    out.push([t[0], m0, m2]);
                    out.push([m0, t[1], m1]);
                    out.push([m2, m1, t[2]]);
                    out.push([m0, m1, m2]);
                }
                2 => {
                    // Rotate so the untouched side is the third one.
                    let r = m.iter().position(|&x| !x).unwrap();
                    let v = [t[(r + 1) % 3], t[(r + 2) % 3], t[r]];
                    let m0 = split(uv, pos, v[0], v[1]);
                    let m1 = split(uv, pos, v[1], v[2]);
                    out.push([m0, v[1], m1]);
                    out.push([v[0], m0, m1]);
                    out.push([v[0], m1, v[2]]);
                }
                1 => {
                    let r = m.iter().position(|&x| x).unwrap();
                    let v = [t[r], t[(r + 1) % 3], t[(r + 2) % 3]];
                    let m0 = split(uv, pos, v[0], v[1]);
                    out.push([v[0], m0, v[2]]);
                    out.push([m0, v[1], v[2]]);
                }
                _ if need[i] => {
                    // Hemmed in by shared edges: the centroid is the only point
                    // this face may add on its own.
                    let p = [
                        (uv[t[0]][0] + uv[t[1]][0] + uv[t[2]][0]) / 3.0,
                        (uv[t[0]][1] + uv[t[1]][1] + uv[t[2]][1]) / 3.0,
                    ];
                    let w = surf.eval(p[0] / mu, p[1] / mv);
                    uv.push(p);
                    pos.push(w);
                    let c = uv.len() - 1;
                    out.push([t[0], t[1], c]);
                    out.push([t[1], t[2], c]);
                    out.push([t[2], t[0], c]);
                }
                _ => out.push(*t),
            }
        }
        *tris = out;
    }
}

// ───────────────────────── polygon triangulation ───────────────────────────

fn cross2(o: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
}

fn area_of(pts: &[[f64; 2]], idx: &[usize]) -> f64 {
    let n = idx.len();
    let mut s = 0.0;
    for i in 0..n {
        let a = pts[idx[i]];
        let b = pts[idx[(i + 1) % n]];
        s += a[0] * b[1] - b[0] * a[1];
    }
    s * 0.5
}

/// Ear clipping, with holes bridged into the outer contour first.
///
/// Returns `None` when the outline cannot be triangulated — a self-crossing
/// boundary, usually from a file the exporter mangled. Reporting that is the
/// point; a fan triangulation would hide it behind a plausible picture.
fn triangulate(pts: &[[f64; 2]], outer: &[usize], holes: &[Vec<usize>]) -> Result<Vec<[usize; 3]>, String> {
    let mut poly: Vec<usize> = outer.to_vec();
    if area_of(pts, &poly) < 0.0 {
        poly.reverse();
    }

    let mut hs: Vec<Vec<usize>> = holes
        .iter()
        .map(|h| {
            let mut h = h.clone();
            if area_of(pts, &h) > 0.0 {
                h.reverse();
            }
            h
        })
        .collect();
    // Bridging works outside-in, so the hole reaching furthest right goes first.
    hs.sort_by(|a, b| {
        let ax = a.iter().map(|&i| pts[i][0]).fold(f64::NEG_INFINITY, f64::max);
        let bx = b.iter().map(|&i| pts[i][0]).fold(f64::NEG_INFINITY, f64::max);
        bx.partial_cmp(&ax).unwrap_or(std::cmp::Ordering::Equal)
    });
    for h in &hs {
        poly = bridge_hole(pts, poly, h)
            .ok_or_else(|| format!("cannot bridge a hole into the outline ({} holes)", hs.len()))?;
    }

    ear_clip(pts, poly).ok_or_else(|| {
        format!(
            "outline crosses itself and could not be split ({} points, {} holes)",
            outer.len(),
            holes.len()
        )
    })
}

/// Ear clipping over one already-bridged ring.
///
/// When no ear is left the ring is **split** on a valid diagonal and each half
/// clipped on its own, which is what carries the awkward cases — a band cut
/// open with a dozen crossing windows bridged into it leaves configurations
/// where every remaining corner is blocked by the far side of a bridge, and
/// giving up there costs the whole face. Only a ring that has no valid diagonal
/// at all comes back `None`.
fn ear_clip(pts: &[[f64; 2]], poly: Vec<usize>) -> Option<Vec<[usize; 3]>> {
    let mut idx = poly;
    if idx.len() < 3 {
        return Some(Vec::new());
    }
    let mut out = Vec::with_capacity(idx.len());
    let mut guard = idx.len() * idx.len() + 16;
    let mut scan = 0usize; // where the last ear was taken from

    // Quality of an ear: area against its longest side, so a fat triangle beats
    // a splinter. Ear clipping is free to take any valid ear, and taking them
    // in blind order produces long blades that reach right across a curved face
    // — geometrically legal in the flat parameter domain, and a sheet cutting
    // through the middle of a tube once mapped back onto it.
    let quality = |a: usize, b: usize, c: usize| -> f64 {
        let (pa, pb, pc) = (pts[a], pts[b], pts[c]);
        let area = cross2(pa, pb, pc).abs();
        let l = |p: [f64; 2], q: [f64; 2]| (p[0] - q[0]).hypot(p[1] - q[1]);
        let longest = l(pa, pb).max(l(pb, pc)).max(l(pc, pa));
        if longest <= 0.0 {
            0.0
        } else {
            area / (longest * longest)
        }
    };

    while idx.len() > 3 {
        let n = idx.len();
        let mut cut: Option<(usize, [usize; 3])> = None;
        for relaxed in [false, true] {
            let mut best = f64::NEG_INFINITY;
            let mut seen = 0usize;
            // Carry on from where the last ear came off rather than restarting
            // at the beginning, which is what makes ear clipping walk round the
            // outline instead of gnawing at one end of it.
            for step in 0..n {
                let i = (scan + step) % n;
                let (a, b, c) = (idx[(i + n - 1) % n], idx[i], idx[(i + 1) % n]);
                if !is_ear(pts, &idx, i, a, b, c, relaxed) {
                    continue;
                }
                let q = quality(a, b, c);
                if q > best {
                    best = q;
                    cut = Some((i, [a, b, c]));
                }
                seen += 1;
                // A handful of candidates is enough to avoid the worst choices
                // without turning every clip into a full scan.
                if seen >= 12 {
                    break;
                }
            }
            if cut.is_some() {
                break;
            }
        }
        let (i, tri) = match cut {
            Some(c) => c,
            None => {
                // Nothing left to bite off: cut the ring in two and clip the
                // halves. Each half is a smaller ring of the same material, so
                // the triangles still tile exactly the region the boundary
                // encloses.
                let (a, b) = split_ring(pts, &idx)?;
                out.extend(ear_clip(pts, a)?);
                out.extend(ear_clip(pts, b)?);
                out.retain(|t| cross2(pts[t[0]], pts[t[1]], pts[t[2]]).abs() > 0.0);
                return Some(out);
            }
        };
        out.push(tri);
        idx.remove(i);
        scan = if i == 0 { 0 } else { i - 1 };
        guard -= 1;
        if guard == 0 {
            return None;
        }
    }
    out.push([idx[0], idx[1], idx[2]]);
    out.retain(|t| cross2(pts[t[0]], pts[t[1]], pts[t[2]]).abs() > 0.0);
    Some(out)
}

/// Cut a ring in two along a diagonal that stays inside it.
///
/// The diagonal must not cross any edge of the ring and must run through the
/// interior rather than round the outside, which is the whole difficulty: on a
/// bridged ring the two copies of a bridge vertex sit on top of each other, so
/// "inside" cannot be judged from the corner alone and the midpoint is tested
/// against the ring itself.
fn split_ring(pts: &[[f64; 2]], idx: &[usize]) -> Option<(Vec<usize>, Vec<usize>)> {
    let n = idx.len();
    let mut best: Option<(f64, usize, usize)> = None;
    for i in 0..n {
        for j in (i + 2)..n {
            if i == 0 && j == n - 1 {
                continue;
            }
            let (pa, pb) = (pts[idx[i]], pts[idx[j]]);
            if pa == pb {
                continue;
            }
            let mid = [(pa[0] + pb[0]) * 0.5, (pa[1] + pb[1]) * 0.5];
            if !point_in_polygon(pts, idx, mid) {
                continue;
            }
            let mut clear = true;
            for k in 0..n {
                let (q1, q2) = (pts[idx[k]], pts[idx[(k + 1) % n]]);
                if segments_cross(pa, pb, q1, q2) {
                    clear = false;
                    break;
                }
            }
            if !clear {
                continue;
            }
            // Prefer a short diagonal: it splits the ring into two pieces that
            // are each still fat, rather than shaving a sliver off one end.
            let d = (pa[0] - pb[0]).hypot(pa[1] - pb[1]);
            if best.map_or(true, |(bd, _, _)| d < bd) {
                best = Some((d, i, j));
            }
        }
    }
    let (_, i, j) = best?;
    let a: Vec<usize> = idx[i..=j].to_vec();
    let mut b: Vec<usize> = idx[j..].to_vec();
    b.extend_from_slice(&idx[..=i]);
    Some((a, b))
}

fn is_ear(pts: &[[f64; 2]], idx: &[usize], at: usize, a: usize, b: usize, c: usize, relaxed: bool) -> bool {
    let (pa, pb, pc) = (pts[a], pts[b], pts[c]);
    let area = cross2(pa, pb, pc);
    if relaxed {
        if area < 0.0 {
            return false;
        }
    } else if area <= 0.0 {
        return false; // reflex or degenerate
    }
    // Vertices are skipped by *index*, not by position in the ring: bridging a
    // hole puts the same vertex into the outline twice, and its other copy sits
    // exactly on a corner of this triangle. Testing it as an interior point
    // rejects every ear next to a bridge, which is how an outline with holes
    // ends up looking self-crossing when it is nothing of the sort.
    for (k, &v) in idx.iter().enumerate() {
        if k == at || v == a || v == b || v == c {
            continue;
        }
        let p = pts[v];
        if p == pa || p == pb || p == pc {
            continue;
        }
        if cross2(pa, pb, p) >= 0.0 && cross2(pb, pc, p) >= 0.0 && cross2(pc, pa, p) >= 0.0 {
            return false;
        }
    }
    true
}

/// Splice a hole into the outer contour with a two-sided bridge.
///
/// The classic construction — cast a ray from the hole's rightmost point and
/// bridge to the first edge it meets — picks the right vertex nearly always,
/// and quietly picks a wrong one often enough to matter on a face with eleven
/// holes in it. A bridge that crosses the outline does not make ear clipping
/// fail; it makes it succeed on a polygon that is no longer simple, and some of
/// the triangles come out inside-out. Those are invisible on screen and show up
/// only as a wrong enclosed volume.
///
/// So every candidate is *checked*: the bridge segment must cross nothing. The
/// ray hit is tried first, then vertices in order of how close they lie to the
/// ray, and the first one that is genuinely visible wins.
fn bridge_hole(pts: &[[f64; 2]], poly: Vec<usize>, hole: &[usize]) -> Option<Vec<usize>> {
    let (hm, _) = hole
        .iter()
        .enumerate()
        .max_by(|a, b| pts[*a.1][0].partial_cmp(&pts[*b.1][0]).unwrap_or(std::cmp::Ordering::Equal))?;
    let m = pts[hole[hm]];
    let n = poly.len();

    // Preferred target: the nearest edge the +x ray meets, taking that edge's
    // outer endpoint (its nearer one when the edge stands straight up).
    let mut first: Option<usize> = None;
    let mut best_x = f64::INFINITY;
    for i in 0..n {
        let a = pts[poly[i]];
        let b = pts[poly[(i + 1) % n]];
        if (a[1] > m[1]) == (b[1] > m[1]) {
            continue;
        }
        let x = a[0] + (m[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
        if x >= m[0] && x < best_x {
            best_x = x;
            first = Some(if a[0] > b[0] {
                i
            } else if a[0] < b[0] {
                (i + 1) % n
            } else if (a[1] - m[1]).abs() <= (b[1] - m[1]).abs() {
                i
            } else {
                (i + 1) % n
            });
        }
    }

    // Everything else, ordered by how far off the ray it sits and then by
    // distance — the same preference the textbook refinement expresses, only
    // applied to the whole outline instead of one triangle.
    let mut rest: Vec<usize> = (0..n).filter(|&i| Some(i) != first).collect();
    rest.sort_by(|&i, &j| {
        let score = |k: usize| {
            let v = pts[poly[k]];
            let dx = v[0] - m[0];
            let ang = if dx > 0.0 { (v[1] - m[1]).abs() / dx } else { f64::INFINITY };
            (ang, (v[0] - m[0]).hypot(v[1] - m[1]))
        };
        let (a1, d1) = score(i);
        let (a2, d2) = score(j);
        a1.partial_cmp(&a2)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(d1.partial_cmp(&d2).unwrap_or(std::cmp::Ordering::Equal))
    });

    let visible = |link: usize| -> bool {
        let p = pts[poly[link]];
        for i in 0..n {
            if i == link || (i + 1) % n == link {
                continue; // edges touching the target meet it legitimately
            }
            if segments_cross(m, p, pts[poly[i]], pts[poly[(i + 1) % n]]) {
                return false;
            }
        }
        for i in 0..hole.len() {
            if i == hm || (i + 1) % hole.len() == hm {
                continue;
            }
            if segments_cross(m, p, pts[hole[i]], pts[hole[(i + 1) % hole.len()]]) {
                return false;
            }
        }
        true
    };

    let link = first
        .filter(|&l| visible(l))
        .or_else(|| rest.iter().copied().find(|&l| visible(l)))
        .or(first)?;

    let mut out = Vec::with_capacity(poly.len() + hole.len() + 2);
    out.extend_from_slice(&poly[..=link]);
    for k in 0..hole.len() {
        out.push(hole[(hm + k) % hole.len()]);
    }
    out.push(hole[hm]);
    out.push(poly[link]);
    out.extend_from_slice(&poly[link + 1..]);
    Some(out)
}

/// True when an outline neither crosses itself nor any of its holes.
fn simple_outline(outer: &[[f64; 2]], holes: &[Vec<[f64; 2]>]) -> bool {
    let n = outer.len();
    for i in 0..n {
        let (a1, a2) = (outer[i], outer[(i + 1) % n]);
        for j in (i + 1)..n {
            if j == i || (j + 1) % n == i || j == (i + 1) % n {
                continue;
            }
            if segments_cross(a1, a2, outer[j], outer[(j + 1) % n]) {
                return false;
            }
        }
        for h in holes {
            for k in 0..h.len() {
                if segments_cross(a1, a2, h[k], h[(k + 1) % h.len()]) {
                    return false;
                }
            }
        }
    }
    true
}

/// Close the small openings a mis-attached loop leaves behind.
///
/// When a face is handed a window that belongs to its neighbour, dropping the
/// window keeps the face — but the *other* face still stops at that window's
/// rim, so the shell has a hole the size of one filament crossing. A hole is
/// fatal to reading the mesh as a solid: a plane through it produces a curve
/// with two loose ends, and a section drawn from that has a line running off
/// across the picture.
///
/// The patch is flat where the true surface is slightly curved, by well under a
/// tenth of a millimetre on a 0.55 mm filament — and it sits *inside* the
/// material, between two filaments that overlap there, so no section and no
/// silhouette ever shows it. What it buys is a mesh that closes.
///
/// Only openings small enough to be one of these is closed. A genuinely open
/// model — a surface that was never a solid — keeps its boundary and is
/// reported rather than quietly lidded.
fn close_small_holes(faces: &mut Vec<FaceMesh>, lo: [f64; 3], hi: [f64; 3]) -> usize {
    if faces.is_empty() {
        return 0;
    }
    const Q: f64 = 1e6;
    let key = |x: f32, y: f32, z: f32| -> (i64, i64, i64) {
        (
            (x as f64 * Q).round() as i64,
            (y as f64 * Q).round() as i64,
            (z as f64 * Q).round() as i64,
        )
    };
    let un = |k: (i64, i64, i64)| [k.0 as f64 / Q, k.1 as f64 / Q, k.2 as f64 / Q];

    // Directed edge counts over welded positions.
    let mut count: HashMap<((i64, i64, i64), (i64, i64, i64)), i32> = HashMap::new();
    let mut solid_of: HashMap<(i64, i64, i64), u32> = HashMap::new();
    for f in faces.iter() {
        for t in f.indices.chunks_exact(3) {
            let k = |i: u32| {
                let a = i as usize * 3;
                key(f.positions[a], f.positions[a + 1], f.positions[a + 2])
            };
            let (ka, kb, kc) = (k(t[0]), k(t[1]), k(t[2]));
            for (x, y) in [(ka, kb), (kb, kc), (kc, ka)] {
                if x == y {
                    continue;
                }
                solid_of.insert(x, f.solid_id);
                *count.entry((x, y)).or_insert(0) += 1;
            }
        }
    }
    // An edge with no partner running the other way is on a boundary.
    let mut next: HashMap<(i64, i64, i64), (i64, i64, i64)> = HashMap::new();
    for (&(a, b), &n) in &count {
        let back = count.get(&(b, a)).copied().unwrap_or(0);
        for _ in 0..(n - back).max(0) {
            next.entry(a).or_insert(b);
        }
    }
    if next.is_empty() {
        return 0;
    }

    let diag =
        ((hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2) + (hi[2] - lo[2]).powi(2)).sqrt();
    let max_span = (diag * 0.05).max(1e-3);

    let mut done: HashSet<(i64, i64, i64)> = HashSet::new();
    let mut capped = 0usize;
    let mut patches: Vec<FaceMesh> = Vec::new();

    let starts: Vec<(i64, i64, i64)> = next.keys().copied().collect();
    for s in starts {
        if done.contains(&s) {
            continue;
        }
        let mut ring: Vec<(i64, i64, i64)> = Vec::new();
        let mut at = s;
        for _ in 0..4096 {
            if !done.insert(at) {
                break;
            }
            ring.push(at);
            match next.get(&at) {
                Some(&n) if n != s => at = n,
                _ => break,
            }
        }
        if ring.len() < 3 {
            continue;
        }
        let pts: Vec<[f64; 3]> = ring.iter().map(|k| un(*k)).collect();
        let mut blo = [f64::MAX; 3];
        let mut bhi = [f64::MIN; 3];
        for p in &pts {
            for a in 0..3 {
                blo[a] = blo[a].min(p[a]);
                bhi[a] = bhi[a].max(p[a]);
            }
        }
        let span = ((bhi[0] - blo[0]).powi(2) + (bhi[1] - blo[1]).powi(2) + (bhi[2] - blo[2]).powi(2))
            .sqrt();
        if span > max_span {
            continue; // a real boundary, not a punched-out window
        }

        // Fan from the centroid: the openings are small and near-elliptical, so
        // a fan is both valid and better shaped than an ear clip would be.
        let mut c = [0.0f64; 3];
        for p in &pts {
            for a in 0..3 {
                c[a] += p[a] / pts.len() as f64;
            }
        }
        // One patch per opening, so each can be turned the right way round on
        // its own — they face every which way and a shared decision would put
        // half of them inside out.
        let mut positions: Vec<f32> = vec![c[0] as f32, c[1] as f32, c[2] as f32];
        let mut normals: Vec<f32> = Vec::new();
        let mut indices: Vec<u32> = Vec::new();
        for p in &pts {
            positions.extend([p[0] as f32, p[1] as f32, p[2] as f32]);
        }
        for i in 0..pts.len() {
            let a = 1 + i as u32;
            let b = 1 + ((i + 1) % pts.len()) as u32;
            indices.extend([0, a, b]);
        }
        // A flat patch: one normal for the lot, from the ring's own plane.
        let mut n = [0.0f64; 3];
        for i in 0..pts.len() {
            let p = pts[i];
            let q = pts[(i + 1) % pts.len()];
            n[0] += (p[1] - c[1]) * (q[2] - c[2]) - (p[2] - c[2]) * (q[1] - c[1]);
            n[1] += (p[2] - c[2]) * (q[0] - c[0]) - (p[0] - c[0]) * (q[2] - c[2]);
            n[2] += (p[0] - c[0]) * (q[1] - c[1]) - (p[1] - c[1]) * (q[0] - c[0]);
        }
        let ln = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt().max(1e-12);
        for _ in 0..=pts.len() {
            normals.extend([
                (n[0] / ln) as f32,
                (n[1] / ln) as f32,
                (n[2] / ln) as f32,
            ]);
        }
        patches.push(FaceMesh {
            face_id: 0,
            solid_id: solid_of.get(&ring[0]).copied().unwrap_or(0),
            positions,
            normals,
            indices,
        });
        capped += 1;
    }

    faces.extend(patches);
    capped
}

/// Make the whole mesh agree about which side is out.
///
/// Each face is triangulated in its surface's own parameter space, and whether
/// that space runs the same way round as the outside of the solid depends on
/// the surface — for a cylinder, on which way its axis was written. Getting it
/// wrong is invisible in a render, which uses the normals, and fatal to
/// anything that reads the mesh as a solid: a plane section takes the side the
/// material is on from the winding, so one face wound backwards puts material
/// on the wrong side of its own wall and the section comes out inside out
/// there.
///
/// Neighbouring faces are made to agree first — two triangles sharing an edge
/// must run along it in opposite directions — and then each connected piece is
/// turned outwards as a whole, by the sign of the volume it encloses. Normals
/// follow the winding afterwards, so the two can no longer disagree.
///
/// Returns the number of triangles that had to be turned round.
fn orient_solid(faces: &mut [FaceMesh]) -> usize {
    if faces.is_empty() {
        return 0;
    }
    // Vertices are matched by position: faces are meshed independently, so the
    // two copies of a shared edge are equal numbers in different buffers.
    const Q: f64 = 1e6; // a millionth of a millimetre
    let mut vid: HashMap<(i64, i64, i64), u32> = HashMap::new();
    // (face, first triangle, triangle count), so a global triangle number can
    // be turned back into somewhere to write.
    let mut span: Vec<(usize, usize, usize)> = Vec::with_capacity(faces.len());
    let mut tri: Vec<[u32; 3]> = Vec::new();
    let mut pos: Vec<[f64; 3]> = Vec::new();
    for (fi, f) in faces.iter().enumerate() {
        let first = tri.len();
        for t in f.indices.chunks_exact(3) {
            let mut w = [0u32; 3];
            for (k, &i) in t.iter().enumerate() {
                let a = i as usize * 3;
                let p = [
                    f.positions[a] as f64,
                    f.positions[a + 1] as f64,
                    f.positions[a + 2] as f64,
                ];
                let key = (
                    (p[0] * Q).round() as i64,
                    (p[1] * Q).round() as i64,
                    (p[2] * Q).round() as i64,
                );
                w[k] = *vid.entry(key).or_insert_with(|| {
                    pos.push(p);
                    (pos.len() - 1) as u32
                });
            }
            tri.push(w);
        }
        span.push((fi, first, tri.len() - first));
    }
    if tri.is_empty() {
        return 0;
    }

    // Every triangle edge, as (low vertex, high vertex, triangle, which way it
    // ran). Sorted rather than hashed: there are three entries per triangle and
    // a hash map of that many keys costs more memory than the mesh does.
    let mut uses: Vec<(u32, u32, u32, bool)> = Vec::with_capacity(tri.len() * 3);
    for (ti, t) in tri.iter().enumerate() {
        for e in 0..3 {
            let (a, b) = (t[e], t[(e + 1) % 3]);
            if a == b {
                continue;
            }
            let fwd = a < b;
            uses.push((a.min(b), a.max(b), ti as u32, fwd));
        }
    }
    uses.sort_unstable();

    // Two triangles sharing an edge agree when they run along it opposite ways.
    let mut adj: Vec<Vec<(u32, bool)>> = vec![Vec::new(); tri.len()];
    let mut i = 0usize;
    while i < uses.len() {
        let mut j = i + 1;
        while j < uses.len() && uses[j].0 == uses[i].0 && uses[j].1 == uses[i].1 {
            j += 1;
        }
        if j - i == 2 {
            let (_, _, ta, da) = uses[i];
            let (_, _, tb, db) = uses[i + 1];
            if ta != tb {
                let same = da == db; // same direction means one is inside out
                adj[ta as usize].push((tb, same));
                adj[tb as usize].push((ta, same));
            }
        }
        i = j;
    }

    // Flood fill: agree with the neighbours, then turn each piece outwards by
    // the sign of the volume it encloses.
    let mut flip = vec![false; tri.len()];
    let mut seen = vec![false; tri.len()];
    let mut part: Vec<u32> = Vec::new();
    let mut queue: Vec<u32> = Vec::new();
    for root in 0..tri.len() {
        if seen[root] {
            continue;
        }
        part.clear();
        queue.clear();
        seen[root] = true;
        queue.push(root as u32);
        part.push(root as u32);
        while let Some(t) = queue.pop() {
            for k in 0..adj[t as usize].len() {
                let (g, must_flip) = adj[t as usize][k];
                if !seen[g as usize] {
                    seen[g as usize] = true;
                    flip[g as usize] = flip[t as usize] ^ must_flip;
                    part.push(g);
                    queue.push(g);
                }
            }
        }
        let mut vol = 0.0f64;
        for &t in &part {
            let w = tri[t as usize];
            let a = pos[w[0] as usize];
            let (b, c) = if flip[t as usize] {
                (pos[w[2] as usize], pos[w[1] as usize])
            } else {
                (pos[w[1] as usize], pos[w[2] as usize])
            };
            vol += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]);
        }
        if vol < 0.0 {
            for &t in &part {
                flip[t as usize] = !flip[t as usize];
            }
        }
    }

    let mut turned = 0usize;
    for (fi, first, count) in span {
        let f = &mut faces[fi];
        for k in 0..count {
            if flip[first + k] {
                f.indices.swap(k * 3 + 1, k * 3 + 2);
                turned += 1;
            }
        }
        // Normals follow the winding, whatever the file said, so that a render
        // and a section can never disagree about which way the wall faces.
        let mut agree = 0.0f64;
        for t in f.indices.chunks_exact(3) {
            let p = |i: u32| {
                let k = i as usize * 3;
                [
                    f.positions[k] as f64,
                    f.positions[k + 1] as f64,
                    f.positions[k + 2] as f64,
                ]
            };
            let (a, b, c) = (p(t[0]), p(t[1]), p(t[2]));
            let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let g = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            for i in t {
                let k = *i as usize * 3;
                agree += g[0] * f.normals[k] as f64
                    + g[1] * f.normals[k + 1] as f64
                    + g[2] * f.normals[k + 2] as f64;
            }
        }
        if agree < 0.0 {
            for n in f.normals.iter_mut() {
                *n = -*n;
            }
        }
    }
    turned
}

/// Does this loop cross the outline rather than sit inside it?
fn loop_crosses(pts: &[[f64; 2]], outer: &[usize], hole: &[usize]) -> bool {
    let n = outer.len();
    let m = hole.len();
    for i in 0..m {
        let (a, b) = (pts[hole[i]], pts[hole[(i + 1) % m]]);
        for j in 0..n {
            if segments_cross(a, b, pts[outer[j]], pts[outer[(j + 1) % n]]) {
                return true;
            }
        }
    }
    false
}

/// Is `q` inside the polygon `poly` (indices into `pts`)?
///
/// Ray crossing to +x, counting edges that straddle the ray in y. Points that
/// land exactly on the outline count as inside, which is what a hole touching
/// the rim wants.
fn point_in_polygon(pts: &[[f64; 2]], poly: &[usize], q: [f64; 2]) -> bool {
    let n = poly.len();
    let mut inside = false;
    for i in 0..n {
        let a = pts[poly[i]];
        let b = pts[poly[(i + 1) % n]];
        if (a[1] > q[1]) != (b[1] > q[1]) {
            let t = (q[1] - a[1]) / (b[1] - a[1]);
            if a[0] + t * (b[0] - a[0]) > q[0] {
                inside = !inside;
            }
        }
    }
    inside
}

/// True when two segments cross properly — sharing an endpoint does not count.
fn segments_cross(p1: [f64; 2], p2: [f64; 2], q1: [f64; 2], q2: [f64; 2]) -> bool {
    let d1 = cross2(q1, q2, p1);
    let d2 = cross2(q1, q2, p2);
    let d3 = cross2(p1, p2, q1);
    let d4 = cross2(p1, p2, q2);
    ((d1 > 0.0) != (d2 > 0.0)) && ((d3 > 0.0) != (d4 > 0.0))
}

// ─────────────────────────────── diagnostics ───────────────────────────────

/// Walk every triangle edge and report the two things that can be wrong with a
/// surface mesh: gaps and inside-out patches.
///
/// In a closed, consistently wound solid every edge is traversed exactly twice,
/// once in each direction. An edge seen only once is a hole. An edge seen twice
/// the *same* way is two triangles that disagree about which side is out. Both
/// are counted rather than inferred from normals, which on a coarse curved
/// triangle can point almost anywhere without anything being wrong.
fn mesh_defects(meshes: &[FaceMesh]) -> (usize, usize) {
    // Positions are quantised to a micron so that vertices which coincide
    // geometrically also coincide as keys.
    let key = |x: f32, y: f32, z: f32| -> (i64, i64, i64) {
        (
            (x as f64 * 1000.0).round() as i64,
            (y as f64 * 1000.0).round() as i64,
            (z as f64 * 1000.0).round() as i64,
        )
    };
    type K = (i64, i64, i64);
    let mut dir: HashMap<(K, K), i32> = HashMap::new();
    for m in meshes {
        for t in m.indices.chunks_exact(3) {
            let v: Vec<K> = t
                .iter()
                .map(|&i| {
                    let o = i as usize * 3;
                    key(m.positions[o], m.positions[o + 1], m.positions[o + 2])
                })
                .collect();
            // A triangle with two coincident corners covers nothing; counting
            // its remaining two edges would pair one of them with itself and
            // read as a disagreement that is not there.
            if v[0] == v[1] || v[1] == v[2] || v[2] == v[0] {
                continue;
            }
            for e in 0..3 {
                *dir.entry((v[e], v[(e + 1) % 3])).or_insert(0) += 1;
            }
        }
    }
    let mut open = 0usize;
    let mut bad = 0usize;
    for (&(a, b), &n) in &dir {
        let back = dir.get(&(b, a)).copied().unwrap_or(0);
        if n + back == 1 {
            open += 1;
        }
        if n > 1 {
            // Counted once per offending pair.
            bad += (n - 1) as usize;
        }
    }
    (open, bad)
}

/// Enclosed volume of a set of face meshes, by the divergence theorem.
///
/// A watertight, consistently wound solid gives its true volume. It is the
/// cheapest check that the import produced the model rather than a plausible
/// pile of surfaces — a stray full torus or a 10 mm phantom tube moves this
/// number immediately.
pub fn mesh_volume(meshes: &[FaceMesh]) -> f64 {
    let mut v = 0.0;
    for m in meshes {
        for t in m.indices.chunks_exact(3) {
            let p: Vec<[f64; 3]> = t
                .iter()
                .map(|&i| {
                    let o = i as usize * 3;
                    [
                        m.positions[o] as f64,
                        m.positions[o + 1] as f64,
                        m.positions[o + 2] as f64,
                    ]
                })
                .collect();
            v += p[0][0] * (p[1][1] * p[2][2] - p[2][1] * p[1][2])
                - p[1][0] * (p[0][1] * p[2][2] - p[2][1] * p[0][2])
                + p[2][0] * (p[0][1] * p[1][2] - p[1][1] * p[0][2]);
        }
    }
    v / 6.0
}
