//! A small, fast, deterministic hasher for the kernel's internal lookup tables.
//!
//! The assembler's vertex weld probes 27 grid cells per point and its shared
//! edge cache is hit once per emitted segment — on a 15 × 15 mm, 12-layer
//! lattice that is millions of lookups of small integer keys.  `std`'s SipHash
//! is DoS-resistant, which nothing here needs, and several times slower.
//!
//! Use it ONLY for maps that are looked up by key and never iterated: the
//! iteration order differs from `std`'s, and a result that depended on it would
//! change.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

/// Multiply–rotate word hasher (FxHash family) with a final avalanche, so the
/// high bits hashbrown uses for its control bytes are well mixed.
#[derive(Default, Clone, Copy)]
pub(crate) struct FastHasher(u64);

const K: u64 = 0x517c_c1b7_2722_0a95;

impl FastHasher {
    #[inline]
    fn add(&mut self, w: u64) {
        self.0 = (self.0.rotate_left(5) ^ w).wrapping_mul(K);
    }
}

impl Hasher for FastHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut chunks = bytes.chunks_exact(8);
        for c in &mut chunks {
            self.add(u64::from_le_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]]));
        }
        let rest = chunks.remainder();
        if !rest.is_empty() {
            let mut b = [0u8; 8];
            b[..rest.len()].copy_from_slice(rest);
            self.add(u64::from_le_bytes(b));
        }
    }
    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.add(i as u64);
    }
    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.add(i as u64);
    }
    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.add(i);
    }
    #[inline]
    fn write_i64(&mut self, i: i64) {
        self.add(i as u64);
    }
    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64);
    }
    #[inline]
    fn finish(&self) -> u64 {
        // murmur3 fmix64
        let mut h = self.0;
        h ^= h >> 33;
        h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
        h ^= h >> 33;
        h = h.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
        h ^ (h >> 33)
    }
}

/// `HashMap` with [`FastHasher`] — lookup tables only, never iterate.
pub(crate) type FastMap<K, V> = HashMap<K, V, BuildHasherDefault<FastHasher>>;

/// `items.iter().map(f).collect()` on all cores, results in input order.
///
/// For pure per-item work on shared read-only data (classifying faces,
/// sectioning them): the output is exactly the serial one.  Items are handed
/// out one at a time, so a few expensive items do not stall a chunk.
pub(crate) fn par_map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get()).min(16);
    if threads <= 1 || items.len() < 16 {
        return items.iter().map(&f).collect();
    }
    let next = std::sync::atomic::AtomicUsize::new(0);
    let mut out: Vec<(usize, R)> = std::thread::scope(|s| {
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                s.spawn(|| {
                    let mut done = Vec::new();
                    loop {
                        let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        if i >= items.len() {
                            break;
                        }
                        done.push((i, f(&items[i])));
                    }
                    done
                })
            })
            .collect();
        workers.into_iter().flat_map(|w| w.join().expect("parallel worker panicked")).collect()
    });
    out.sort_by_key(|(i, _)| *i);
    out.into_iter().map(|(_, r)| r).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn par_map_keeps_order() {
        let v: Vec<u64> = (0..1000).collect();
        let r = par_map(&v, |x| x * x);
        assert_eq!(r, v.iter().map(|x| x * x).collect::<Vec<_>>());
    }

    #[test]
    fn grid_keys_spread_and_round_trip() {
        let mut m: FastMap<(i64, i64, i64), usize> = FastMap::default();
        let mut n = 0;
        for x in -20..20 {
            for y in -20..20 {
                for z in -5..5 {
                    m.insert((x, y, z), n);
                    n += 1;
                }
            }
        }
        assert_eq!(m.len(), n);
        let mut n2 = 0;
        for x in -20..20 {
            for y in -20..20 {
                for z in -5..5 {
                    assert_eq!(m.get(&(x, y, z)), Some(&n2));
                    n2 += 1;
                }
            }
        }
    }
}
