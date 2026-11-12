//! Filament fixtures captured from `brep::pipeline` (CADCORE_DUMP_FILAMENTS)
//! — real generator output, used by the universal-fuse tests.
#![cfg(test)]
#![allow(clippy::excessive_precision)]

use super::filament::{ElbowSeg, LegSeg};
use super::scaffold::Filament;
use cadcore_geom::{CylSurf, TorusSurf};
use cadcore_math::{Point3, UnitVec3, Vec3};

pub(crate) mod gen {
    use super::*;
    /// 3 mm × 2 layers (the `fuse_then_cut` pipeline fixture).
    pub(crate) fn pair_3mm() -> Vec<Filament> {
        let mut out = Vec::new();
        {
            let legs = vec![
                LegSeg { surf: CylSurf::new(Point3::new(0.299999982, 0.275000006, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.000000000, 1.000000000, 0.000000000)).unwrap(), 0.275000000), length: 1.958232987 },
                LegSeg { surf: CylSurf::new(Point3::new(0.775144693, 2.724186258, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.999464628, 0.032717846, 0.000000000)).unwrap(), 0.275000000), length: 0.216918685 },
                LegSeg { surf: CylSurf::new(Point3::new(1.498972047, 2.256401702, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.032717966, -0.999464624, 0.000000000)).unwrap(), 0.275000000), length: 1.508086781 },
                LegSeg { surf: CylSurf::new(Point3::new(2.023194798, 0.274240662, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.999464594, -0.032718875, 0.000000000)).unwrap(), 0.275000000), length: 0.155998835 },
                LegSeg { surf: CylSurf::new(Point3::new(2.686386654, 0.756692489, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.006916129, 0.999976083, 0.000000000)).unwrap(), 0.275000000), length: 1.968354492 },
            ];
            let elbows = vec![
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.791216231, 2.233232993, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.491216249, 0.275000000), theta_lo: 1.570796326795, theta_hi: 3.108868967138 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(1.008018783, 2.240330105, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.491216249, 0.275000000), theta_lo: 3.108868967138, theta_hi: 4.679665174728 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(2.039266842, 0.765193911, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.491216249, 0.275000000), theta_lo: 1.603520132452, theta_hi: 3.108867937564 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(2.195182154, 0.760089804, 0.349999994), UnitVec3::try_from_vec(Vec3::new(-0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.491216249, 0.275000000), theta_lo: 3.108867937564, theta_hi: 4.705472796212 },
            ];
            out.push(Filament { legs, elbows });
        }
        {
            let legs = vec![
                LegSeg { surf: CylSurf::new(Point3::new(2.724999905, 2.700000048, 0.699999988), UnitVec3::try_from_vec(Vec3::new(-1.000000000, 0.000000000, 0.000000000)).unwrap(), 0.275000000), length: 1.958232640 },
                LegSeg { surf: CylSurf::new(Point3::new(0.275813654, 2.224855492, 0.699999988), UnitVec3::try_from_vec(Vec3::new(-0.032718875, -0.999464594, 0.000000000)).unwrap(), 0.275000000), length: 0.216918003 },
                LegSeg { surf: CylSurf::new(Point3::new(0.743598343, 1.501027948, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.999464624, -0.032717966, 0.000000000)).unwrap(), 0.275000000), length: 1.508086254 },
                LegSeg { surf: CylSurf::new(Point3::new(2.725759223, 0.976804358, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.032717846, -0.999464628, 0.000000000)).unwrap(), 0.275000000), length: 0.155997726 },
                LegSeg { surf: CylSurf::new(Point3::new(2.243306914, 0.313613738, 0.699999988), UnitVec3::try_from_vec(Vec3::new(-0.999976082, -0.006916315, 0.000000000)).unwrap(), 0.275000000), length: 1.968354017 },
            ];
            let elbows = vec![
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.766767264, 2.208783437, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.491216611, 0.275000000), theta_lo: 0.000000000000, theta_hi: 1.538071610769 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.759669951, 1.991981573, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.491216611, 0.275000000), theta_lo: 1.538071610769, theta_hi: 3.108868847933 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(2.234805596, 0.960732809, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.491216611, 0.275000000), theta_lo: -3.108868847933, theta_hi: -1.603520013247 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(2.239909505, 0.804818599, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.000000000, -0.000000000, -1.000000000)).unwrap(), 0.491216611, 0.275000000), theta_lo: -1.603520013247, theta_hi: -0.006916369748 },
            ];
            out.push(Filament { legs, elbows });
        }
        out
    }
    /// 6 mm × 4 layers orthogonal grid, pitch 1.2, D 0.55, bend 0.5.
    pub(crate) fn grid_6mm_4l() -> Vec<Filament> {
        let mut out = Vec::new();
        {
            let legs = vec![
                LegSeg { surf: CylSurf::new(Point3::new(0.275000006, 0.275000006, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.000000000, 1.000000000, 0.000000000)).unwrap(), 0.275000006), length: 5.048440686 },
                LegSeg { surf: CylSurf::new(Point3::new(0.660298687, 5.721556637, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.999464794, 0.032712775, 0.000000000)).unwrap(), 0.275000006), length: 0.290061099 },
                LegSeg { surf: CylSurf::new(Point3::new(1.361350842, 5.345962827, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.032720256, -0.999464549, 0.000000000)).unwrap(), 0.275000006), length: 4.689844655 },
                LegSeg { surf: CylSurf::new(Point3::new(1.899886545, 0.273546926, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.999464579, -0.032719329, 0.000000000)).unwrap(), 0.275000006), length: 0.140390555 },
                LegSeg { surf: CylSurf::new(Point3::new(2.451350842, 0.654035893, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.032720195, 0.999464551, 0.000000000)).unwrap(), 0.275000006), length: 4.689845515 },
                LegSeg { surf: CylSurf::new(Point3::new(2.989888876, 5.726452772, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.999464792, 0.032712836, 0.000000000)).unwrap(), 0.275000006), length: 0.140390729 },
                LegSeg { surf: CylSurf::new(Point3::new(3.541350789, 5.345962838, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.032720256, -0.999464549, 0.000000000)).unwrap(), 0.275000006), length: 4.689844633 },
                LegSeg { surf: CylSurf::new(Point3::new(4.079886861, 0.273546947, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.999464610, -0.032718401, 0.000000000)).unwrap(), 0.275000006), length: 0.140391268 },
                LegSeg { surf: CylSurf::new(Point3::new(4.631351602, 0.654038989, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.032712775, 0.999464794, 0.000000000)).unwrap(), 0.275000006), length: 4.689837206 },
                LegSeg { surf: CylSurf::new(Point3::new(5.169854662, 5.726451657, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.999464792, 0.032712836, 0.000000000)).unwrap(), 0.275000006), length: 0.128428955 },
                LegSeg { surf: CylSurf::new(Point3::new(5.709572637, 5.333751736, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.003049605, -0.999995350, 0.000000000)).unwrap(), 0.275000006), length: 5.058775164 },
            ];
            let elbows = vec![
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.673329139, 5.323440692, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398329133, 0.275000006), theta_lo: 1.570796326795, theta_hi: 3.108874040966 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.963234995, 5.332929396, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398329133, 0.275000006), theta_lo: 3.108874040966, theta_hi: 4.679662883008 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(1.912919607, 0.671662785, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.398329133, 0.275000006), theta_lo: 1.603522424172, theta_hi: 3.108867484213 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(2.053234994, 0.667069300, 0.349999994), UnitVec3::try_from_vec(Vec3::new(-0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.398329133, 0.275000006), theta_lo: 3.108867484213, theta_hi: 4.679662943769 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(3.002919351, 5.328336828, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398329133, 0.275000006), theta_lo: 1.603522363411, theta_hi: 3.108873980218 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(3.143234942, 5.332929407, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398329133, 0.275000006), theta_lo: 3.108873980218, theta_hi: 4.679662883008 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(4.092919554, 0.671662818, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.398329133, 0.275000006), theta_lo: 1.603522424172, theta_hi: 3.108868412214 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(4.233235658, 0.667069440, 0.349999994), UnitVec3::try_from_vec(Vec3::new(-0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.398329133, 0.275000006), theta_lo: 3.108868412214, theta_hi: 4.679670367761 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(5.182885138, 5.328335714, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398329133, 0.275000006), theta_lo: 1.603514939419, theta_hi: 3.108873980218 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(5.311245356, 5.332536989, 0.349999994), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398329133, 0.275000006), theta_lo: 3.108873980218, theta_hi: 4.709339370290 },
            ];
            out.push(Filament { legs, elbows });
        }
        {
            let legs = vec![
                LegSeg { surf: CylSurf::new(Point3::new(5.724999905, 0.275000006, 0.699999988), UnitVec3::try_from_vec(Vec3::new(-1.000000000, 0.000000000, 0.000000000)).unwrap(), 0.275000006), length: 5.048442656 },
                LegSeg { surf: CylSurf::new(Point3::new(0.278444015, 0.660293941, 0.699999988), UnitVec3::try_from_vec(Vec3::new(-0.032718280, 0.999464614, 0.000000000)).unwrap(), 0.275000006), length: 0.290066322 },
                LegSeg { surf: CylSurf::new(Point3::new(0.654033457, 1.361350731, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.999464551, 0.032720195, 0.000000000)).unwrap(), 0.275000006), length: 4.689850425 },
                LegSeg { surf: CylSurf::new(Point3::new(5.726452691, 1.899886373, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.032712836, 0.999464792, 0.000000000)).unwrap(), 0.275000006), length: 0.140395960 },
                LegSeg { surf: CylSurf::new(Point3::new(5.345965481, 2.451350788, 0.699999988), UnitVec3::try_from_vec(Vec3::new(-0.999464549, 0.032720256, 0.000000000)).unwrap(), 0.275000006), length: 4.689849724 },
                LegSeg { surf: CylSurf::new(Point3::new(0.273547051, 2.989884571, 0.699999988), UnitVec3::try_from_vec(Vec3::new(-0.032718158, 0.999464618, 0.000000000)).unwrap(), 0.275000006), length: 0.140395555 },
                LegSeg { surf: CylSurf::new(Point3::new(0.654033482, 3.541350679, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.999464551, 0.032720195, 0.000000000)).unwrap(), 0.275000006), length: 4.689850408 },
                LegSeg { surf: CylSurf::new(Point3::new(5.726452699, 4.079886320, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.032712836, 0.999464792, 0.000000000)).unwrap(), 0.275000006), length: 0.140396668 },
                LegSeg { surf: CylSurf::new(Point3::new(5.345962556, 4.631351541, 0.699999988), UnitVec3::try_from_vec(Vec3::new(-0.999464792, 0.032712836, 0.000000000)).unwrap(), 0.275000006), length: 4.689841645 },
                LegSeg { surf: CylSurf::new(Point3::new(0.273548121, 5.169849891, 0.699999988), UnitVec3::try_from_vec(Vec3::new(-0.032719329, 0.999464579, 0.000000000)).unwrap(), 0.275000006), length: 0.128433807 },
                LegSeg { surf: CylSurf::new(Point3::new(0.666244328, 5.709572547, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.999995350, 0.003049621, 0.000000000)).unwrap(), 0.275000006), length: 5.058779100 },
            ];
            let elbows = vec![
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.676557249, 0.673326498, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398326492, 0.275000006), theta_lo: 0.000000000000, theta_hi: 1.538072206934 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.667066778, 0.963237522, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398326492, 0.275000006), theta_lo: 1.538072206934, theta_hi: 3.108866616974 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(5.328339386, 1.912916762, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.398326492, 0.275000006), theta_lo: -3.108866616974, theta_hi: -1.603515000166 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(5.332932136, 2.053237580, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.000000000, -0.000000000, 1.000000000)).unwrap(), 0.398326492, 0.275000006), theta_lo: -1.603515000166, theta_hi: -0.032726097377 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.671660286, 3.002917080, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398326492, 0.275000006), theta_lo: 0.032726097377, theta_hi: 1.538072328449 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.667066802, 3.143237471, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398326492, 0.275000006), theta_lo: 1.538072328449, theta_hi: 3.108866616974 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(5.328339394, 4.092916710, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.398326492, 0.275000006), theta_lo: -3.108866616974, theta_hi: -1.603515000166 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(5.332932167, 4.233238236, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.000000000, -0.000000000, 1.000000000)).unwrap(), 0.398326492, 0.275000006), theta_lo: -1.603515000166, theta_hi: -0.032718673371 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.671661341, 5.182882867, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398326492, 0.275000006), theta_lo: 0.032718673371, theta_hi: 1.538071157418 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.667459073, 5.311247907, 0.699999988), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398326492, 0.275000006), theta_lo: 1.538071157418, theta_hi: 3.138543028124 },
            ];
            out.push(Filament { legs, elbows });
        }
        {
            let legs = vec![
                LegSeg { surf: CylSurf::new(Point3::new(5.724999905, 5.724999905, 1.049999952), UnitVec3::try_from_vec(Vec3::new(0.000000000, -1.000000000, 0.000000000)).unwrap(), 0.275000006), length: 5.048451092 },
                LegSeg { surf: CylSurf::new(Point3::new(5.339714925, 0.278444414, 1.049999952), UnitVec3::try_from_vec(Vec3::new(-0.999464579, -0.032719329, 0.000000000)).unwrap(), 0.275000006), length: 0.290084786 },
                LegSeg { surf: CylSurf::new(Point3::new(4.638648285, 0.654027381, 1.049999952), UnitVec3::try_from_vec(Vec3::new(-0.032712715, 0.999464796, 0.000000000)).unwrap(), 0.275000006), length: 4.689859516 },
                LegSeg { surf: CylSurf::new(Point3::new(4.100155864, 5.726451297, 1.049999952), UnitVec3::try_from_vec(Vec3::new(-0.999464792, 0.032712836, 0.000000000)).unwrap(), 0.275000006), length: 0.140447522 },
                LegSeg { surf: CylSurf::new(Point3::new(3.548649024, 5.345974306, 1.049999952), UnitVec3::try_from_vec(Vec3::new(-0.032720256, -0.999464549, 0.000000000)).unwrap(), 0.275000006), length: 4.689866847 },
                LegSeg { surf: CylSurf::new(Point3::new(3.010123310, 0.273547298, 1.049999952), UnitVec3::try_from_vec(Vec3::new(-0.999464610, -0.032718401, 0.000000000)).unwrap(), 0.275000006), length: 0.140412409 },
                LegSeg { surf: CylSurf::new(Point3::new(2.458649371, 0.654024632, 1.049999952), UnitVec3::try_from_vec(Vec3::new(-0.032720135, 0.999464553, 0.000000000)).unwrap(), 0.275000006), length: 4.689867479 },
                LegSeg { surf: CylSurf::new(Point3::new(1.920121988, 5.726452411, 1.049999952), UnitVec3::try_from_vec(Vec3::new(-0.999464792, 0.032712836, 0.000000000)).unwrap(), 0.275000006), length: 0.140413751 },
                LegSeg { surf: CylSurf::new(Point3::new(1.368648852, 5.345972838, 1.049999952), UnitVec3::try_from_vec(Vec3::new(-0.032716546, -0.999464671, 0.000000000)).unwrap(), 0.275000006), length: 4.689862764 },
                LegSeg { surf: CylSurf::new(Point3::new(0.830140623, 0.273547863, 1.049999952), UnitVec3::try_from_vec(Vec3::new(-0.999464610, -0.032718401, 0.000000000)).unwrap(), 0.275000006), length: 0.128431280 },
                LegSeg { surf: CylSurf::new(Point3::new(0.290429974, 0.666235287, 1.049999952), UnitVec3::try_from_vec(Vec3::new(-0.003050131, 0.999995348, 0.000000000)).unwrap(), 0.275000006), length: 5.058788150 },
            ];
            let elbows = vec![
                ElbowSeg { surf: TorusSurf::new(Point3::new(5.326682239, 0.676548813, 1.049999952), UnitVec3::try_from_vec(Vec3::new(0.000000000, -0.000000000, -1.000000000)).unwrap(), 0.398317666, 0.275000006), theta_lo: -1.570796326795, theta_hi: -0.032725169376 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(5.036752770, 0.667057433, 1.049999952), UnitVec3::try_from_vec(Vec3::new(-0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398317666, 0.275000006), theta_lo: -0.032725169376, theta_hi: 1.538077774918 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(4.087125764, 5.328346813, 1.049999952), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.398317666, 0.275000006), theta_lo: -1.538077774918, theta_hi: -0.032718673371 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(3.946753410, 5.332941250, 1.049999952), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.398317666, 0.275000006), theta_lo: -0.032718673371, theta_hi: 1.538070229418 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(2.997090993, 0.671651709, 1.049999952), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398317666, 0.275000006), theta_lo: -1.538070229418, theta_hi: -0.032724241376 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(2.856753759, 0.667057639, 1.049999952), UnitVec3::try_from_vec(Vec3::new(-0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398317666, 0.275000006), theta_lo: -0.032724241376, theta_hi: 1.538070350940 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(1.907091888, 5.328347928, 1.049999952), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.398317666, 0.275000006), theta_lo: -1.538070350940, theta_hi: -0.032718673371 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(1.766753287, 5.332941260, 1.049999952), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.398317666, 0.275000006), theta_lo: -0.032718673371, theta_hi: 1.538073941420 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.817108306, 0.671652274, 1.049999952), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398317666, 0.275000006), theta_lo: -1.538073941420, theta_hi: -0.032724241376 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.688745787, 0.667450208, 1.049999952), UnitVec3::try_from_vec(Vec3::new(-0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398317666, 0.275000006), theta_lo: -0.032724241376, theta_hi: 1.567746190775 },
            ];
            out.push(Filament { legs, elbows });
        }
        {
            let legs = vec![
                LegSeg { surf: CylSurf::new(Point3::new(0.275000006, 5.724999905, 1.399999976), UnitVec3::try_from_vec(Vec3::new(1.000000000, 0.000000000, 0.000000000)).unwrap(), 0.275000006), length: 5.048449683 },
                LegSeg { surf: CylSurf::new(Point3::new(5.721556334, 5.339710247, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.032712836, -0.999464792, 0.000000000)).unwrap(), 0.275000006), length: 0.290080503 },
                LegSeg { surf: CylSurf::new(Point3::new(5.345969216, 4.638648181, 1.399999976), UnitVec3::try_from_vec(Vec3::new(-0.999464792, -0.032712836, 0.000000000)).unwrap(), 0.275000006), length: 4.689854509 },
                LegSeg { surf: CylSurf::new(Point3::new(0.273548354, 4.100155483, 1.399999976), UnitVec3::try_from_vec(Vec3::new(-0.032718401, -0.999464610, 0.000000000)).unwrap(), 0.275000006), length: 0.140442753 },
                LegSeg { surf: CylSurf::new(Point3::new(0.654026737, 3.548648986, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.999464549, -0.032720256, 0.000000000)).unwrap(), 0.275000006), length: 4.689863440 },
                LegSeg { surf: CylSurf::new(Point3::new(5.726452504, 3.010119076, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.032712836, -0.999464792, 0.000000000)).unwrap(), 0.275000006), length: 0.140408228 },
                LegSeg { surf: CylSurf::new(Point3::new(5.345972137, 2.458649277, 1.399999976), UnitVec3::try_from_vec(Vec3::new(-0.999464549, -0.032720256, 0.000000000)).unwrap(), 0.275000006), length: 4.689862679 },
                LegSeg { surf: CylSurf::new(Point3::new(0.273547214, 1.920121930, 1.399999976), UnitVec3::try_from_vec(Vec3::new(-0.032719207, -0.999464583, 0.000000000)).unwrap(), 0.275000006), length: 0.140408981 },
                LegSeg { surf: CylSurf::new(Point3::new(0.654028057, 1.368648821, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.999464671, -0.032716546, 0.000000000)).unwrap(), 0.275000006), length: 4.689859504 },
                LegSeg { surf: CylSurf::new(Point3::new(5.726451937, 0.830136367, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.032712775, -0.999464794, 0.000000000)).unwrap(), 0.275000006), length: 0.128427052 },
                LegSeg { surf: CylSurf::new(Point3::new(5.333761421, 0.290429927, 1.399999976), UnitVec3::try_from_vec(Vec3::new(-0.999995348, -0.003050124, 0.000000000)).unwrap(), 0.275000006), length: 5.058784858 },
            ];
            let elbows = vec![
                ElbowSeg { surf: TorusSurf::new(Point3::new(5.323449689, 5.326680076, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398319829, 0.275000006), theta_lo: -3.141592653590, theta_hi: -1.603515000166 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(5.332939045, 5.036754826, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.000000000, -0.000000000, -1.000000000)).unwrap(), 0.398319829, 0.275000006), theta_lo: -1.603515000166, theta_hi: -0.032718673371 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.671654926, 4.087123095, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.398319829, 0.275000006), theta_lo: 0.032718673371, theta_hi: 1.538072085419 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.667059864, 3.946755534, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.398319829, 0.275000006), theta_lo: 1.538072085419, theta_hi: 3.108866556213 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(5.328345859, 2.997088905, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398319829, 0.275000006), theta_lo: -3.108866556213, theta_hi: -1.603515000166 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(5.332939010, 2.856755825, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.000000000, -0.000000000, -1.000000000)).unwrap(), 0.398319829, 0.275000006), theta_lo: -1.603515000166, theta_hi: -0.032726097377 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.671653776, 1.907089221, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.398319829, 0.275000006), theta_lo: 0.032726097377, theta_hi: 1.538071278937 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(0.667059706, 1.766755417, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, 1.000000000)).unwrap(), 0.398319829, 0.275000006), theta_lo: 1.538071278937, theta_hi: 3.108870268215 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(5.328345291, 0.817106220, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.000000000, 0.000000000, -1.000000000)).unwrap(), 0.398319829, 0.275000006), theta_lo: -3.108870268215, theta_hi: -1.603514939419 },
                ElbowSeg { surf: TorusSurf::new(Point3::new(5.332546496, 0.688747903, 1.399999976), UnitVec3::try_from_vec(Vec3::new(0.000000000, -0.000000000, -1.000000000)).unwrap(), 0.398319829, 0.275000006), theta_lo: -1.603514939419, theta_hi: -0.003050128696 },
            ];
            out.push(Filament { legs, elbows });
        }
        out
    }
}

/// Load filaments from a `CADCORE_DUMP_FILAMENTS` log (the `[fildump]` lines of
/// `brep::pipeline`), so any generator output can be replayed against the
/// kernel without pasting it into the source.
pub(crate) fn from_dump(path: &str) -> Vec<Filament> {
    let text = std::fs::read_to_string(path).expect("read filament dump");
    let nums = |l: &str| -> Vec<f64> {
        let l = l.replace("Point3", "").replace("UnitVec3", "").replace("Vec3", "");
        let mut out = Vec::new();
        let mut cur = String::new();
        for ch in l.chars().chain(std::iter::once(' ')) {
            if ch.is_ascii_digit() || ch == '.' || ch == '-' || ch == 'e' && !cur.is_empty() {
                cur.push(ch);
            } else if !cur.is_empty() {
                if let Ok(v) = cur.parse::<f64>() {
                    out.push(v);
                }
                cur.clear();
            }
        }
        out
    };
    let mut out: Vec<Filament> = Vec::new();
    let mut legs: Vec<LegSeg> = Vec::new();
    let mut elbows: Vec<ElbowSeg> = Vec::new();
    let mut open = false;
    let mut flush = |legs: &mut Vec<LegSeg>, elbows: &mut Vec<ElbowSeg>, out: &mut Vec<Filament>| {
        if !legs.is_empty() {
            out.push(Filament { legs: std::mem::take(legs), elbows: std::mem::take(elbows) });
        }
    };
    for line in text.lines() {
        let Some(l) = line.split("[fildump]").nth(1) else { continue };
        if l.contains("// filament") {
            flush(&mut legs, &mut elbows, &mut out);
            open = true;
        } else if l.contains("LegSeg {") && open {
            let v = nums(l);
            assert_eq!(v.len(), 8, "leg line: {l}");
            legs.push(LegSeg {
                surf: CylSurf::new(
                    Point3::new(v[0], v[1], v[2]),
                    UnitVec3::try_from_vec(Vec3::new(v[3], v[4], v[5])).unwrap(),
                    v[6],
                ),
                length: v[7],
            });
        } else if l.contains("ElbowSeg {") && open {
            let v = nums(l);
            assert_eq!(v.len(), 10, "elbow line: {l}");
            elbows.push(ElbowSeg {
                surf: TorusSurf::new(
                    Point3::new(v[0], v[1], v[2]),
                    UnitVec3::try_from_vec(Vec3::new(v[3], v[4], v[5])).unwrap(),
                    v[6],
                    v[7],
                ),
                theta_lo: v[8],
                theta_hi: v[9],
            });
        }
    }
    flush(&mut legs, &mut elbows, &mut out);
    out
}
