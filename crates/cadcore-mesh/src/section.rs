//! The outline a plane cuts out of a solid, as closed curves the figure can draw.
//!
//! Cutting a triangle mesh with a plane is easy and gives the wrong picture. A
//! fused scaffold does not arrive as one solid: the exporter writes it as a
//! stack of Z-band shells, and consecutive bands **overlap** — that is what
//! makes the print continuous. Chain the raw cut segments and every filament
//! crossing comes out as two rectangles drawn straight through each other,
//! which is exactly what a section must never show. The material there is one
//! body; the drawing has to say so.
//!
//! So the cut is treated as a plane region rather than as a set of lines. Each
//! segment is oriented from its triangle's own winding, so that material lies
//! on its left; the segments are split wherever they cross; and a piece is kept
//! only where the material stops — where the winding number on its right-hand
//! side is zero. Pieces buried inside another filament fail that test and drop
//! out, and what is left is the boundary of the union: pore outlines, and the
//! outside of the block. Nothing is approximated on a grid and nothing is
//! smoothed, so the result is exactly as accurate as the mesh it came from.
//!
//! The same rule handles the other two cases without a special path: one clean
//! solid (nothing overlaps, so nothing is dropped) and a solid touching itself
//! (the shared wall is dropped from both sides).

/// Knobs for a section. The defaults suit a 0.55 mm filament drawn at page size.
#[derive(Clone, Copy, Debug)]
pub struct SectionOptions {
    /// Distance below which two cut points are the same point, in millimetres.
    /// Only ever matches points that *should* coincide: adjacent triangles
    /// interpolate the same shared edge from opposite ends, which agrees to the
    /// last bit or two and to nothing coarser.
    pub weld: f64,
    /// Douglas–Peucker tolerance, in millimetres. The cut of a tessellated
    /// cylinder repeats a straight line a thousand times; dropping the points
    /// that lie on it changes no pixel and makes the union tractable. Zero
    /// keeps every point.
    pub simplify: f64,
    /// Take the union of everything the plane cuts. Off draws the raw cut,
    /// crossings and all, which is only useful for seeing what the union did.
    pub union: bool,
}

impl Default for SectionOptions {
    fn default() -> Self {
        Self {
            weld: 1e-7,
            simplify: 2e-3,
            union: true,
        }
    }
}

/// What the plane cut, and how much work it took.
#[derive(Clone, Debug, Default)]
pub struct SectionOutline {
    /// Closed curves in the plane's own two coordinates, in millimetres. The
    /// two are the model axes the section axis leaves out, in ascending order:
    /// a Z section is (x, y), a Y section is (x, z), an X section is (y, z).
    pub loops: Vec<Vec<[f64; 2]>>,
    /// Whether each curve in `loops` came back to its own start. A curve that
    /// did not is where the body has a gap, and it must be drawn as an open
    /// line rather than closed with an invented one.
    pub closed: Vec<bool>,
    /// Triangle-plane segments the cut produced, before any of the rest.
    pub segments: usize,
    /// Chains that did not close. Zero on a watertight mesh; anything else is
    /// a hole in the model, and the chain is still returned so it can be seen.
    pub open_chains: usize,
    /// Pieces the union threw away as interior — the crossings.
    pub interior_dropped: usize,
}

/// The two model axes a section along `axis` keeps, in ascending order.
pub fn plane_axes(axis: usize) -> (usize, usize) {
    match axis {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    }
}

/// Cut `positions`/`indices` with the plane `axis = position` and return the
/// boundary of the material as closed curves.
pub fn section_outline(
    positions: &[f32],
    indices: &[u32],
    axis: usize,
    position: f64,
    o: &SectionOptions,
) -> SectionOutline {
    let (a1, a2) = plane_axes(axis);
    let mut segs = cut_triangles(positions, indices, axis, position, a1, a2);
    let mut out = SectionOutline {
        segments: segs.len(),
        ..Default::default()
    };
    if segs.is_empty() {
        return out;
    }

    // Which way round the cut runs is decided once for the whole of it. Two of
    // the three plane frames are left-handed against the model's, which would
    // put the material on the right; the total signed area says which, and it
    // takes the whole cut to say so rather than one sample that a single
    // inside-out triangle could mislead.
    Welder::new(o.weld).snap_all(&mut segs);
    if segs
        .iter()
        .map(|(a, b)| a[0] * b[1] - b[0] * a[1])
        .sum::<f64>()
        < 0.0
    {
        for s in &mut segs {
            *s = (s.1, s.0);
        }
    }

    // The union runs on the cut segments themselves, before anything is threaded
    // into curves.
    //
    // Threading first looked natural and is a trap. Where the outline touches
    // itself — a doubled wall, two filaments meeting at a point — several
    // segments leave the same place and any threading has to guess; guessing
    // wrong leaves a curve in pieces. None of that troubles the segments: every
    // vertex still has as many arriving as leaving, so the winding number is
    // exactly defined whatever a threader would have made of them. Threading is
    // left to the end, where the worst it can cost is an extra `M` in the path.
    if !o.union {
        let mut raw = segs.clone();
        let mut open = 0usize;
        let chains = chain(&mut raw, o.weld, &mut open);
        out.open_chains = open;
        out.loops = chains.iter().map(|c| c.pts.clone()).collect();
        out.closed = chains.iter().map(|c| c.closed).collect();
        return out;
    }

    let (kept, dropped) = union_boundary(&segs, o.weld);
    out.interior_dropped = dropped;
    let mut pieces = kept;
    let mut chains = chain(&mut pieces, o.weld, &mut out.open_chains);

    if o.simplify > 0.0 {
        for c in &mut chains {
            c.pts = if c.closed {
                simplify_closed(&c.pts, o.simplify)
            } else {
                dp(&c.pts, o.simplify)
            };
        }
    }
    for c in &mut chains {
        if c.closed {
            drop_collinear(&mut c.pts, o.weld);
        }
    }
    chains.retain(|c| c.pts.len() >= if c.closed { 3 } else { 2 });
    out.loops = chains.iter().map(|c| c.pts.clone()).collect();
    out.closed = chains.iter().map(|c| c.closed).collect();
    out
}

// ─────────────────────────────── plane × mesh ────────────────────────────────

/// One oriented cut segment, in plane coordinates.
type Seg = ([f64; 2], [f64; 2]);

/// A run of segments read head to tail, and whether it came back to its start.
///
/// The distinction matters more than it looks. A closed run bounds an area, so
/// it can be filled and it contributes its wrap-around edge to the winding
/// count. An open one is where the body has a gap; joining its ends with a
/// straight line would invent a wall across the picture and — worse — put that
/// invented wall into the winding count, which then reports material on the
/// wrong side for everything in its shadow.
struct Chain {
    pts: Vec<[f64; 2]>,
    closed: bool,
}

/// Where the plane meets each triangle, oriented so the solid is on the left.
///
/// The direction is `m × n`: `m` is the plane normal and `n` the triangle's own
/// outward normal, so the segment runs along the wall with the wall's outside
/// on the right. Taking it from the winding rather than from a containment test
/// is what lets the union work on a mesh that overlaps itself.
fn cut_triangles(
    p: &[f32],
    idx: &[u32],
    axis: usize,
    position: f64,
    a1: usize,
    a2: usize,
) -> Vec<Seg> {
    let mut segs: Vec<Seg> = Vec::new();
    let at = |i: u32, k: usize| p[i as usize * 3 + k] as f64;
    let d = |i: u32| at(i, axis) - position;

    // Which side of the plane a corner is on, with "exactly on it" counted as
    // below. Some rule is needed and this one is the only kind that works:
    // *consistent*. Treating a corner on the plane as a third state gives a
    // triangle with an edge in the plane two crossing points at the same place,
    // and both triangles sharing that edge emit it, so the cut comes back with
    // doubled and zero-length pieces and nothing chains. Sections through a
    // printed lattice land on such planes constantly — the midplane between two
    // layers is where the tessellation puts a whole row of vertices — so this is
    // the normal case, not a corner case.
    let above = |d: f64| d > 0.0;

    for t in idx.chunks_exact(3) {
        let (ia, ib, ic) = (t[0], t[1], t[2]);
        let (da, db, dc) = (d(ia), d(ib), d(ic));
        let (sa, sb, sc) = (above(da), above(db), above(dc));
        if sa == sb && sb == sc {
            continue;
        }

        let mut hit: [[f64; 3]; 2] = [[0.0; 3]; 2];
        let mut n = 0usize;
        let edge = |i: u32, j: u32, di: f64, dj: f64, hit: &mut [[f64; 3]; 2], n: &mut usize| {
            if above(di) == above(dj) || *n >= 2 {
                return;
            }
            let s = di / (di - dj);
            hit[*n] = [
                at(i, 0) + (at(j, 0) - at(i, 0)) * s,
                at(i, 1) + (at(j, 1) - at(i, 1)) * s,
                at(i, 2) + (at(j, 2) - at(i, 2)) * s,
            ];
            *n += 1;
        };
        edge(ia, ib, da, db, &mut hit, &mut n);
        edge(ib, ic, db, dc, &mut hit, &mut n);
        edge(ic, ia, dc, da, &mut hit, &mut n);
        if n < 2 {
            continue;
        }
        let (p0, p1) = (hit[0], hit[1]);
        if p0 == p1 {
            continue;
        }

        // Outward normal from the winding, and the direction the cut should run.
        let (pa, pb, pc) = (
            [at(ia, 0), at(ia, 1), at(ia, 2)],
            [at(ib, 0), at(ib, 1), at(ib, 2)],
            [at(ic, 0), at(ic, 1), at(ic, 2)],
        );
        let u = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
        let v = [pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]];
        let nor = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        // m x n with m the axis unit vector.
        let mut want = [0.0f64; 3];
        want[(axis + 1) % 3] = -nor[(axis + 2) % 3];
        want[(axis + 2) % 3] = nor[(axis + 1) % 3];
        let dir = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let along = want[0] * dir[0] + want[1] * dir[1] + want[2] * dir[2];

        let q0 = [p0[a1], p0[a2]];
        let q1 = [p1[a1], p1[a2]];
        if along >= 0.0 {
            segs.push((q0, q1));
        } else {
            segs.push((q1, q0));
        }
    }
    segs
}

// ──────────────────────────────── chaining ───────────────────────────────────

fn qkey(p: [f64; 2], weld: f64) -> (i64, i64) {
    let q = 1.0 / weld.max(1e-12);
    ((p[0] * q).round() as i64, (p[1] * q).round() as i64)
}

/// Snap points that should be the same point to one shared value.
///
/// Two things make near-duplicates: adjacent triangles interpolate a shared
/// edge from opposite ends, and a crossing is solved once along each of the two
/// edges that meet there. Both agree to the last bit or two, and both must come
/// out **identical**, because everything downstream chains by exact equality.
/// Rounding to a grid alone is not enough — two values either side of a cell
/// boundary round apart — so each point also looks in the neighbouring cells
/// for a representative already claimed.
struct Welder {
    weld: f64,
    seen: std::collections::HashMap<(i64, i64), [f64; 2]>,
}

impl Welder {
    fn new(weld: f64) -> Self {
        Self {
            weld,
            seen: std::collections::HashMap::new(),
        }
    }

    fn snap(&mut self, p: [f64; 2]) -> [f64; 2] {
        let k = qkey(p, self.weld);
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(q) = self.seen.get(&(k.0 + dx, k.1 + dy)) {
                    if (q[0] - p[0]).abs() <= self.weld && (q[1] - p[1]).abs() <= self.weld {
                        return *q;
                    }
                }
            }
        }
        self.seen.insert(k, p);
        p
    }

    fn snap_all(&mut self, segs: &mut [Seg]) {
        for s in segs.iter_mut() {
            s.0 = self.snap(s.0);
            s.1 = self.snap(s.1);
        }
    }
}

/// Drop points that lie on the straight line between their neighbours.
///
/// Splitting an edge at a T-junction leaves a vertex in the middle of a
/// straight run; so does welding two collinear pieces back together. They are
/// harmless but they show up as extra nodes in an exported SVG, and a figure
/// meant to be edited in Illustrator should not carry them.
fn drop_collinear(l: &mut Vec<[f64; 2]>, tol: f64) {
    if l.len() < 4 {
        return;
    }
    let mut out: Vec<[f64; 2]> = Vec::with_capacity(l.len());
    let n = l.len();
    for i in 0..n {
        let a = l[(i + n - 1) % n];
        let b = l[i];
        let c = l[(i + 1) % n];
        let cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        let base = (c[0] - a[0]).hypot(c[1] - a[1]);
        if base > 0.0 && (cross / base).abs() <= tol {
            continue;
        }
        out.push(b);
    }
    if out.len() >= 3 {
        *l = out;
    }
}

/// Thread oriented segments head-to-tail into closed curves.
///
/// Direction matters: a segment is only ever followed by one that *starts*
/// where it ends. Matching either way round, as an unoriented chainer must,
/// jumps between two filaments wherever their cuts happen to touch, and welds
/// the whole section into one meaningless curve.
fn chain(segs: &mut Vec<Seg>, weld: f64, open: &mut usize) -> Vec<Chain> {
    use std::collections::HashMap;
    let mut start: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (i, s) in segs.iter().enumerate() {
        start.entry(qkey(s.0, weld)).or_default().push(i);
    }
    let mut used = vec![false; segs.len()];
    let mut out: Vec<Chain> = Vec::new();

    let dir = |s: &Seg| {
        let (dx, dy) = (s.1[0] - s.0[0], s.1[1] - s.0[1]);
        let l = dx.hypot(dy).max(1e-15);
        [dx / l, dy / l]
    };

    for i in 0..segs.len() {
        if used[i] {
            continue;
        }
        used[i] = true;
        let head = segs[i].0;
        let mut pts = vec![segs[i].0, segs[i].1];
        let mut came = dir(&segs[i]);
        let mut closed = false;
        loop {
            let tail = *pts.last().unwrap();
            if qkey(tail, weld) == qkey(head, weld) && pts.len() > 2 {
                pts.pop(); // the closing point is implicit
                closed = true;
                break;
            }
            // Where the outline touches itself — a doubled wall, or two solids
            // meeting at a point — several pieces leave the same place, and
            // taking whichever came first tears one curve into several. The
            // continuation that turns least is the one that was the same curve.
            let mut best: Option<(f64, usize)> = None;
            for &j in start.get(&qkey(tail, weld)).map_or(&[][..], |v| &v[..]) {
                if used[j] {
                    continue;
                }
                let d = dir(&segs[j]);
                let turn = came[0] * d[0] + came[1] * d[1];
                if best.map_or(true, |(b, _)| turn > b) {
                    best = Some((turn, j));
                }
            }
            match best {
                Some((_, j)) => {
                    used[j] = true;
                    came = dir(&segs[j]);
                    pts.push(segs[j].1);
                }
                None => break,
            }
        }
        if !closed {
            *open += 1;
        }
        if pts.len() >= 2 {
            out.push(Chain { pts, closed });
        }
    }
    out
}

#[cfg_attr(not(test), allow(dead_code))]
fn signed_area(l: &[[f64; 2]]) -> f64 {
    let n = l.len();
    let mut s = 0.0;
    for i in 0..n {
        let a = l[i];
        let b = l[(i + 1) % n];
        s += a[0] * b[1] - b[0] * a[1];
    }
    s * 0.5
}

// ─────────────────────────────── simplifying ─────────────────────────────────

/// Douglas–Peucker on a closed curve, anchored at its two furthest-apart points
/// so the result does not depend on where the curve happened to start.
fn simplify_closed(l: &[[f64; 2]], tol: f64) -> Vec<[f64; 2]> {
    let n = l.len();
    if n < 8 {
        return l.to_vec();
    }
    // Anchor: the point furthest from the first, and the first itself.
    let mut far = 0usize;
    let mut best = -1.0;
    for (i, p) in l.iter().enumerate() {
        let d = (p[0] - l[0][0]).hypot(p[1] - l[0][1]);
        if d > best {
            best = d;
            far = i;
        }
    }
    let mut a: Vec<[f64; 2]> = l[..=far].to_vec();
    let mut b: Vec<[f64; 2]> = l[far..].to_vec();
    b.push(l[0]);
    let mut out = dp(&a, tol);
    out.pop();
    let mut second = dp(&b, tol);
    second.pop();
    out.append(&mut second);
    a.clear();
    out
}

fn dp(pts: &[[f64; 2]], tol: f64) -> Vec<[f64; 2]> {
    if pts.len() < 3 {
        return pts.to_vec();
    }
    let (a, b) = (pts[0], pts[pts.len() - 1]);
    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
    let len = dx.hypot(dy);
    let mut worst = 0.0;
    let mut at = 0usize;
    for (i, p) in pts.iter().enumerate().take(pts.len() - 1).skip(1) {
        let d = if len < 1e-15 {
            (p[0] - a[0]).hypot(p[1] - a[1])
        } else {
            ((p[0] - a[0]) * dy - (p[1] - a[1]) * dx).abs() / len
        };
        if d > worst {
            worst = d;
            at = i;
        }
    }
    if worst <= tol {
        return vec![a, b];
    }
    let mut left = dp(&pts[..=at], tol);
    let right = dp(&pts[at..], tol);
    left.pop();
    left.extend(right);
    left
}

// ───────────────────────────────── the union ─────────────────────────────────

/// Split every edge where it crosses another, then keep only the pieces with
/// nothing on their right — the boundary of the material.
fn union_boundary(segs: &[Seg], weld: f64) -> (Vec<Seg>, usize) {
    if segs.is_empty() {
        return (Vec::new(), 0);
    }
    let grid = Grid::build(segs);
    let mut stamp = vec![u32::MAX; segs.len()];
    let mut pieces: Vec<Seg> = Vec::with_capacity(segs.len() + segs.len() / 4);
    let mut ts: Vec<f64> = Vec::new();
    for (i, e) in segs.iter().enumerate() {
        ts.clear();
        grid.each_near(e, &mut stamp, i as u32, |j| {
            if j == i {
                return;
            }
            match cross_param(*e, segs[j], weld) {
                Cross::At(t) => ts.push(t),
                // Walls lying on top of one another are cut at each other's
                // ends, so the overlapping stretch comes out as one piece on
                // both and can be recognised as the same wall stated twice.
                Cross::Along(t0, t1) => {
                    ts.push(t0);
                    ts.push(t1);
                }
                Cross::None => {}
            }
        });
        ts.retain(|t| *t > 1e-12 && *t < 1.0 - 1e-12);
        if ts.is_empty() {
            pieces.push(*e);
            continue;
        }
        ts.push(0.0);
        ts.push(1.0);
        ts.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for w in ts.windows(2) {
            if w[1] - w[0] < 1e-12 {
                continue;
            }
            pieces.push((lerp(*e, w[0]), lerp(*e, w[1])));
        }
    }

    // The two segments meeting at a crossing each solve for it separately and
    // agree only to the last bit. Welding here is what lets the kept pieces
    // thread back into closed curves by exact equality.
    Welder::new(weld).snap_all(&mut pieces);

    // Where does the material stop? Where the winding number on the right-hand
    // side of a piece is zero — nothing on that side — and the left is inside.
    // Pieces buried in another filament have material on both sides and drop
    // out, which is exactly the crossings.
    let wind = Grid::build(&pieces);
    // How far to step off a piece to ask what is on that side.
    //
    // Small: the step has to stay inside whatever thin sliver of material or
    // air happens to be there, and a section through a lattice is full of
    // places where two walls all but touch. A tenth of a micron on a
    // twenty-millimetre part clears the arithmetic and crosses nothing.
    let eps = (grid.diagonal() * 1e-7).clamp(1e-9, 1e-4);
    let mut kept: Vec<Seg> = Vec::with_capacity(pieces.len());
    let mut dropped = 0usize;
    for s in &pieces {
        let mid = [(s.0[0] + s.1[0]) * 0.5, (s.0[1] + s.1[1]) * 0.5];
        let (dx, dy) = (s.1[0] - s.0[0], s.1[1] - s.0[1]);
        let len = dx.hypot(dy);
        if len < 1e-14 {
            continue;
        }
        let right = [mid[0] + dy / len * eps, mid[1] - dx / len * eps];
        if wind.winding(&pieces, right) == 0 {
            kept.push(*s);
        } else {
            dropped += 1;
        }
    }

    // Coincident walls. Two solids meeting face to face state the same wall
    // twice facing opposite ways — an interior wall, and neither copy belongs
    // on the union. Two overlapping along a wall state it twice facing the same
    // way — one boundary stated twice, and exactly one copy belongs. Both only
    // become visible now, because the collinear split above is what makes the
    // duplicates identical.
    let removed = dedupe_walls(&mut kept, weld);
    (kept, dropped + removed)
}

fn lerp(e: Seg, t: f64) -> [f64; 2] {
    [
        e.0[0] + (e.1[0] - e.0[0]) * t,
        e.0[1] + (e.1[1] - e.0[1]) * t,
    ]
}

/// How one edge meets another, in the first edge's own parameter.
enum Cross {
    None,
    /// A proper crossing at this parameter.
    At(f64),
    /// The two lie on the same line and overlap between these parameters.
    Along(f64, f64),
}

fn cross_param(a: Seg, b: Seg, weld: f64) -> Cross {
    let r = [a.1[0] - a.0[0], a.1[1] - a.0[1]];
    let s = [b.1[0] - b.0[0], b.1[1] - b.0[1]];
    let denom = r[0] * s[1] - r[1] * s[0];
    let len = r[0].hypot(r[1]);
    if len < 1e-15 {
        return Cross::None;
    }
    if (denom / len).abs() <= weld {
        // Parallel. Collinear only if b's ends sit on a's line.
        let off = |q: [f64; 2]| ((q[0] - a.0[0]) * r[1] - (q[1] - a.0[1]) * r[0]) / len;
        if off(b.0).abs() > weld || off(b.1).abs() > weld {
            return Cross::None;
        }
        let along = |q: [f64; 2]| ((q[0] - a.0[0]) * r[0] + (q[1] - a.0[1]) * r[1]) / (len * len);
        let (t0, t1) = {
            let (u, v) = (along(b.0), along(b.1));
            (u.min(v), u.max(v))
        };
        if t1 <= 0.0 || t0 >= 1.0 {
            return Cross::None;
        }
        return Cross::Along(t0.clamp(0.0, 1.0), t1.clamp(0.0, 1.0));
    }
    let qp = [b.0[0] - a.0[0], b.0[1] - a.0[1]];
    let t = (qp[0] * s[1] - qp[1] * s[0]) / denom;
    let u = (qp[0] * r[1] - qp[1] * r[0]) / denom;
    if !(0.0..=1.0).contains(&u) {
        return Cross::None;
    }
    Cross::At(t)
}

fn dedupe_walls(segs: &mut Vec<Seg>, weld: f64) -> usize {
    use std::collections::HashMap;
    type Key = ((i64, i64), (i64, i64));
    let key = |s: &Seg| -> Key { (qkey(s.0, weld), qkey(s.1, weld)) };
    let mut same: HashMap<Key, Vec<usize>> = HashMap::new();
    for (i, s) in segs.iter().enumerate() {
        same.entry(key(s)).or_default().push(i);
    }
    let mut drop = vec![false; segs.len()];
    let mut removed = 0usize;
    for (i, s) in segs.iter().enumerate() {
        if drop[i] {
            continue;
        }
        let back = (qkey(s.1, weld), qkey(s.0, weld));
        if let Some(j) = same.get(&back).and_then(|l| l.iter().copied().find(|&j| !drop[j])) {
            drop[i] = true;
            drop[j] = true;
            removed += 2;
            continue;
        }
        // Any further copies pointing the same way are the same wall restated.
        if let Some(l) = same.get(&key(s)) {
            for &j in l {
                if j != i && !drop[j] {
                    drop[j] = true;
                    removed += 1;
                }
            }
        }
    }
    let mut k = 0usize;
    segs.retain(|_| {
        let d = drop[k];
        k += 1;
        !d
    });
    removed
}

// ───────────────────────── a row index over the segments ─────────────────────

/// Segments bucketed into horizontal bands, which is all the acceleration a
/// section needs: crossings are found band by band, and the winding ray runs
/// along a row so it only ever meets one band's worth of candidates.
struct Grid {
    y0: f64,
    dy: f64,
    rows: Vec<Vec<usize>>,
    bbox: [f64; 4],
}

impl Grid {
    fn build(segs: &[Seg]) -> Self {
        let mut bbox = [f64::MAX, f64::MAX, f64::MIN, f64::MIN];
        for s in segs {
            for p in [s.0, s.1] {
                bbox[0] = bbox[0].min(p[0]);
                bbox[1] = bbox[1].min(p[1]);
                bbox[2] = bbox[2].max(p[0]);
                bbox[3] = bbox[3].max(p[1]);
            }
        }
        let n = segs.len().max(1);
        let rows = ((n as f64).sqrt() as usize).clamp(1, 4096);
        let span = (bbox[3] - bbox[1]).max(1e-9);
        let dy = span / rows as f64;
        let mut g = Grid {
            y0: bbox[1],
            dy,
            rows: vec![Vec::new(); rows],
            bbox,
        };
        for (i, s) in segs.iter().enumerate() {
            let (lo, hi) = g.row_span(s);
            for r in lo..=hi {
                g.rows[r].push(i);
            }
        }
        g
    }

    fn diagonal(&self) -> f64 {
        (self.bbox[2] - self.bbox[0]).hypot(self.bbox[3] - self.bbox[1])
    }

    fn row_of(&self, y: f64) -> usize {
        let r = ((y - self.y0) / self.dy).floor();
        (r.max(0.0) as usize).min(self.rows.len() - 1)
    }

    fn row_span(&self, s: &Seg) -> (usize, usize) {
        let a = self.row_of(s.0[1].min(s.1[1]));
        let b = self.row_of(s.0[1].max(s.1[1]));
        (a.min(b), a.max(b))
    }

    fn each_near(&self, s: &Seg, stamp: &mut [u32], tick: u32, mut f: impl FnMut(usize)) {
        let (lo, hi) = self.row_span(s);
        for r in lo..=hi {
            for &j in &self.rows[r] {
                if stamp[j] == tick {
                    continue;
                }
                stamp[j] = tick;
                f(j);
            }
        }
    }

    /// Winding number of `q`, from the segments crossing a ray to +x.
    ///
    /// The half-open rule — a segment counts when its lower endpoint is at or
    /// below the ray and its upper endpoint strictly above — makes a vertex on
    /// the ray count exactly once however many edges meet there, so no jitter
    /// or perturbation is needed.
    fn winding(&self, segs: &[Seg], q: [f64; 2]) -> i32 {
        let mut w = 0i32;
        for &i in &self.rows[self.row_of(q[1])] {
            let (a, b) = segs[i];
            if a[1] <= q[1] {
                if b[1] > q[1] && side(a, b, q) > 0.0 {
                    w += 1;
                }
            } else if b[1] <= q[1] && side(a, b, q) < 0.0 {
                w -= 1;
            }
        }
        w
    }
}

/// Which side of the line a→b the point q is on.
fn side(a: [f64; 2], b: [f64; 2], q: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (q[1] - a[1]) - (b[1] - a[1]) * (q[0] - a[0])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An axis-aligned box as a triangle soup with outward winding.
    fn cube(lo: [f64; 3], hi: [f64; 3]) -> (Vec<f32>, Vec<u32>) {
        let c = [
            [lo[0], lo[1], lo[2]],
            [hi[0], lo[1], lo[2]],
            [hi[0], hi[1], lo[2]],
            [lo[0], hi[1], lo[2]],
            [lo[0], lo[1], hi[2]],
            [hi[0], lo[1], hi[2]],
            [hi[0], hi[1], hi[2]],
            [lo[0], hi[1], hi[2]],
        ];
        let quads: [[usize; 4]; 6] = [
            [0, 3, 2, 1], // -z
            [4, 5, 6, 7], // +z
            [0, 1, 5, 4], // -y
            [3, 7, 6, 2], // +y
            [0, 4, 7, 3], // -x
            [1, 2, 6, 5], // +x
        ];
        let mut p: Vec<f32> = Vec::new();
        for v in &c {
            p.extend(v.iter().map(|x| *x as f32));
        }
        let mut idx: Vec<u32> = Vec::new();
        for q in &quads {
            idx.extend([q[0] as u32, q[1] as u32, q[2] as u32]);
            idx.extend([q[0] as u32, q[2] as u32, q[3] as u32]);
        }
        (p, idx)
    }

    fn area(o: &SectionOutline) -> f64 {
        o.loops.iter().map(|l| signed_area(l)).sum()
    }

    #[test]
    fn one_box_cuts_to_its_own_rectangle() {
        let (p, i) = cube([0.0, 0.0, 0.0], [2.0, 3.0, 1.0]);
        for axis in 0..3 {
            let mid = [1.0, 1.5, 0.5][axis];
            let o = section_outline(&p, &i, axis, mid, &SectionOptions::default());
            assert_eq!(o.loops.len(), 1, "axis {axis}");
            assert_eq!(o.open_chains, 0, "axis {axis}");
            let want = [3.0 * 1.0, 2.0 * 1.0, 2.0 * 3.0][axis];
            assert!(
                (area(&o) - want).abs() < 1e-9,
                "axis {axis}: {} want {want}",
                area(&o)
            );
        }
    }

    #[test]
    fn two_crossing_bars_draw_one_outline_not_two() {
        // A plus sign: the union has twelve corners and the area of the two
        // bars less their shared middle. Chained raw, this is two rectangles
        // crossing — which is the bug the union exists to fix.
        let (p1, i1) = cube([0.0, 1.0, 0.0], [6.0, 2.0, 1.0]);
        let (p2, i2) = cube([2.0, -1.0, 0.0], [3.0, 4.0, 1.0]);
        let mut p = p1.clone();
        p.extend_from_slice(&p2);
        let mut i = i1.clone();
        i.extend(i2.iter().map(|k| k + (p1.len() / 3) as u32));

        let o = section_outline(&p, &i, 2, 0.5, &SectionOptions::default());
        assert_eq!(o.open_chains, 0);
        assert_eq!(o.loops.len(), 1, "the plus sign is one curve");
        assert_eq!(o.loops[0].len(), 12, "a plus sign has twelve corners");
        let want = 6.0 * 1.0 + 1.0 * 5.0 - 1.0 * 1.0;
        assert!((area(&o) - want).abs() < 1e-9, "{} want {want}", area(&o));
        assert!(o.interior_dropped > 0, "the crossing pieces must be dropped");
    }

    #[test]
    fn a_hole_stays_a_hole() {
        // A square frame: outer 6x6, inner 2x2 void. The void's boundary is
        // material on one side and nothing on the other, so it survives — and
        // it must come back wound the other way.
        let (p, i) = frame();
        let o = section_outline(&p, &i, 2, 0.5, &SectionOptions::default());
        assert_eq!(o.loops.len(), 2);
        let mut a: Vec<f64> = o.loops.iter().map(|l| signed_area(l)).collect();
        a.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert!((a[0] + 4.0).abs() < 1e-9, "inner loop {} want -4", a[0]);
        assert!((a[1] - 36.0).abs() < 1e-9, "outer loop {} want 36", a[1]);
    }

    #[test]
    fn overlapping_boxes_lose_the_shared_wall() {
        let (p1, i1) = cube([0.0, 0.0, 0.0], [2.0, 2.0, 1.0]);
        let (p2, i2) = cube([1.0, 0.0, 0.0], [3.0, 2.0, 1.0]);
        let mut p = p1.clone();
        p.extend_from_slice(&p2);
        let mut i = i1.clone();
        i.extend(i2.iter().map(|k| k + (p1.len() / 3) as u32));
        let o = section_outline(&p, &i, 2, 0.5, &SectionOptions::default());
        assert_eq!(o.loops.len(), 1);
        assert!((area(&o) - 6.0).abs() < 1e-9, "{} want 6", area(&o));
    }

    /// 6x6 square with a 2x2 hole through it, as a closed triangle soup.
    fn frame() -> (Vec<f32>, Vec<u32>) {
        // Built from four boxes round the void, which overlap at the corners —
        // a harder input than a hand-made frame and closer to the real thing.
        let parts = [
            ([0.0, 0.0, 0.0], [6.0, 2.0, 1.0]),
            ([0.0, 4.0, 0.0], [6.0, 6.0, 1.0]),
            ([0.0, 0.0, 0.0], [2.0, 6.0, 1.0]),
            ([4.0, 0.0, 0.0], [6.0, 6.0, 1.0]),
        ];
        let mut p: Vec<f32> = Vec::new();
        let mut i: Vec<u32> = Vec::new();
        for (lo, hi) in parts {
            let (pp, ii) = cube(lo, hi);
            let base = (p.len() / 3) as u32;
            p.extend_from_slice(&pp);
            i.extend(ii.iter().map(|k| k + base));
        }
        (p, i)
    }
}
