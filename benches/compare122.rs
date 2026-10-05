//! Field, hash, and representation costs over
//! GF(2^122) = GF(2^61)\[u\]/(u^2 + u + 1).
//!
//!   nix develop -c cargo bench --bench compare122
//!
//! The dense-constant and subfield-constant curve fixtures come from
//! `tests/common/kats122.rs`. This suite measures point representations
//! complementary to `group.rs`, and direct hashing to the λ-affine and
//! unscaled addends, by try-and-increment and by Pornin's map, which
//! avoids a preparation inversion.
//! `group.rs`'s `binary-lambda.122` hash rows instead return ordinary affine
//! points and share the `binary.122` hash.

mod common;
use common::{each, whole};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use ephemeral_ecmh::curve::binary::{Curve, Pornin, sum_batch, unscaled};
use ephemeral_ecmh::curve::binary122::mul_u2;
use ephemeral_ecmh::curve::h2c::map1;
use ephemeral_ecmh::ecmh::{TAG_ITEM, digest_batch};
use ephemeral_ecmh::field::gf2_122::{self, Gf, Gf61, from_u128};
use ephemeral_ecmh::hash::{Salted, halves};
use std::hint::black_box;

const N: usize = 1024;
/// Independent chains for the throughput-bound products.
const ACCS: usize = 8;

/// N pseudo-random field elements, from digests as hash-to-curve sees them.
fn xs() -> Vec<Gf> {
    let h = Salted::new(b"bench-field", &[0; 32]);
    (0..N as u32)
        .map(|i| from_u128(halves(&h.digest(&[], i))[0]))
        .collect()
}

fn field(c: &mut Criterion) {
    let gx = xs();
    let mut g = c.benchmark_group("field");
    // Latency: one chain of dependent steps. Throughput: ACCS chains. The
    // step count is opaque: mul_u2 is linear of order 3, and a known count
    // lets the compiler fold a chain of it to count mod 3 steps.
    macro_rules! chain {
        ($name:expr, $v:expr, |$a:ident| $step:expr) => {
            g.throughput(Throughput::Elements(64));
            g.bench_function(concat!($name, " latency (dependent chain)"), |bn| {
                bn.iter(|| (0..black_box(64)).fold(black_box($v), |$a, _| $step))
            });
            g.bench_function(concat!($name, " throughput (8 chains)"), |bn| {
                bn.iter(|| {
                    // opaque as a whole: lanes known to be equal could share one chain
                    let mut acc = black_box([$v; ACCS]);
                    for _ in 0..black_box(64 / ACCS) {
                        for lane in acc.iter_mut() {
                            let $a = *lane;
                            *lane = $step;
                        }
                    }
                    acc
                })
            });
        };
    }
    let w = gx[1];
    chain!("gf2_122/mul", gx[0], |a| a * w);
    // squarings chain in a sqrt or an inversion's addition chain
    chain!("gf2_122/square", gx[0], |a| a.square());
    // m_beta on the GLS curves, whose constant is in GF(2^61): 2 base
    // products, where mul is 3
    let beta: Gf61 = gx[1].parts().0;
    chain!("gf2_122/mul_base (by a GF(2^61) constant)", gx[0], |a| a
        .mul_base(beta));
    // a^2 = 1 + a = u^2 for a = u: mul_a2 in the extended and unscaled
    // adds, mul_1a in prepare and normalize; a word swap and an XOR
    chain!("gf2_122/mul_u2 (by a^2 = 1 + a)", gx[0], |a| mul_u2(a));
    each(&mut g, "gf2_122/invert", &gx, |v| v.invert());
    each(&mut g, "gf2_122/sqrt", &gx, |v| v.sqrt());
    // what the halftrace does for odd-degree fields: solve z^2 + z = c
    each(
        &mut g,
        "gf2_122/qsolve (z^2 + z = c, 2 base halftraces)",
        &gx,
        |v| v.qsolve(),
    );
    each(&mut g, "gf2_122/trace", &gx, |v| v.trace());
    each(&mut g, "gf2_122/normalize (to_u128)", &gx, |&v| {
        gf2_122::to_u128(v)
    });
    whole(&mut g, "gf2_122/batch invert (product tree)", &gx, |v| {
        let mut w = v.to_vec();
        ephemeral_ecmh::field::batch::invert(&mut w);
        w
    });
    g.finish();
}

/// One family's rows: `name` is its id prefix, `lambda` that of the same
/// curves under λ accumulators.
fn curve<M: Pornin<F = Gf>>(
    c: &mut Criterion,
    name: &str,
    lambda: &str,
    unscaled: &str,
    fixture: common::families::Fixture<Curve<M>>,
) {
    let bin = fixture.group;
    let salt = Salted::new(TAG_ITEM, &fixture.seed);
    let items = common::items(&[], N);
    let refs: Vec<&[u8]> = items.iter().map(|x| x.as_slice()).collect();
    let hb = bin.hash_to_curve_batch(&salt, &refs);
    let eb: Vec<_> = hb.iter().map(|p| bin.from_affine(p)).collect();

    // Pornin's map, as `compare.rs` times it over GF(2^127); and the λ
    // families' own hash: (x, λ) out, no lift.
    let mut g = c.benchmark_group("hash_to_curve");
    each(&mut g, &format!("{name}/pornin map x1"), &refs, |m| {
        map1(&bin, &salt, m)
    });
    whole(
        &mut g,
        &format!("{name}/pornin map x1, batched"),
        &refs,
        |ms| bin.hash_to_curve_map1_batch(&salt, ms),
    );
    each(
        &mut g,
        &format!("{lambda}/try-and-increment to (x, λ)"),
        &refs,
        |m| bin.hash_to_lambda(&salt, m),
    );
    whole(
        &mut g,
        &format!("{lambda}/try-and-increment to (x, λ), batched"),
        &refs,
        |ms| bin.hash_to_lambda_batch(&salt, ms),
    );
    // The map straight to a prepared addend, as `compare.rs` and
    // `compare109.rs` time it: hash and prepare in one row, no lift.
    let halves_of =
        |ms: &[&[u8]]| -> Vec<u128> { ms.iter().map(|m| halves(&salt.digest(m, 0))[0]).collect() };
    each(
        &mut g,
        &format!("{lambda}/pornin map x1 to (x, λ)"),
        &refs,
        |m| bin.map_to_lambda(halves(&salt.digest(m, 0))[0]),
    );
    whole(
        &mut g,
        &format!("{lambda}/pornin map x1 to (x, λ), batched"),
        &refs,
        |ms| bin.map_to_lambda_batch(&halves_of(ms)),
    );
    // the w codec's addend is the same (x, λ): its own rows, since the
    // report keys hash-to-addend recipes by family
    let wcodec = lambda.replace("-lambda", "-w");
    each(
        &mut g,
        &format!("{wcodec}/pornin map x1 to (x, λ)"),
        &refs,
        |m| bin.map_to_lambda(halves(&salt.digest(m, 0))[0]),
    );
    whole(
        &mut g,
        &format!("{wcodec}/pornin map x1 to (x, λ), batched"),
        &refs,
        |ms| bin.map_to_lambda_batch(&halves_of(ms)),
    );
    let u = unscaled::Curve::new(bin);
    each(
        &mut g,
        &format!("{unscaled}/pornin map x1 to (u, v)"),
        &refs,
        |m| u.map_to_addend(halves(&salt.digest(m, 0))[0]),
    );
    whole(
        &mut g,
        &format!("{unscaled}/pornin map x1 to (u, v), batched"),
        &refs,
        |ms| u.map_to_addend_batch(&halves_of(ms)),
    );
    g.finish();

    let mut g = c.benchmark_group("add");
    g.throughput(Throughput::Elements(N as u64));
    // the affine law needs no curve constant: once, on the dense family
    if name == "gf2_122" {
        g.bench_function(format!("{name}/batch affine (tree sum)"), |bn| {
            bn.iter(|| sum_batch(black_box(&hb)))
        });
    }
    // throughput-bound: ACCS accumulators round robin, like updating
    // distinct RIBLT cells
    g.bench_function(
        format!("{name}/extended += affine, {ACCS} accumulators"),
        |bn| {
            bn.iter(|| {
                let mut acc = [bin.neutral(); ACCS];
                for ch in hb.chunks(ACCS) {
                    for (a, p) in acc.iter_mut().zip(ch) {
                        *a = bin.add_affine(a, p);
                    }
                }
                acc
            })
        },
    );
    whole(&mut g, &format!("{name}/normalize, batched"), &eb, |ps| {
        bin.normalize_batch(ps)
    });
    g.finish();

    let mut g = c.benchmark_group("digest");
    g.throughput(Throughput::Elements(N as u64));
    g.bench_function(format!("{name}/batch"), |bn| {
        bn.iter(|| digest_batch(bin, &fixture.seed, &refs))
    });
    g.finish();
}

fn curves(c: &mut Criterion) {
    curve(
        c,
        "gf2_122",
        "gf2_122-lambda",
        "gf2_122-u",
        common::families::binary122(),
    );
    curve(
        c,
        "gf2_122-gls",
        "gf2_122-gls-lambda",
        "gf2_122-gls-u",
        common::families::binary122_gls(),
    );
    // the affine negation, the same on both; the addends' is group.neg in
    // benches/group.rs
    let k = common::families::binary122();
    let bin = k.group;
    let items = common::items(&[], N);
    let refs: Vec<&[u8]> = items.iter().map(|x| x.as_slice()).collect();
    let hb = bin.hash_to_curve_batch(&Salted::new(TAG_ITEM, &k.seed), &refs);
    let mut g = c.benchmark_group("negate");
    each(&mut g, "gf2_122/affine (y += x)", &hb, |p| p.neg());
    g.finish();
}

criterion_group! {
    name = benches;
    config = common::config();
    targets = field, curves
}
criterion_main!(benches);
