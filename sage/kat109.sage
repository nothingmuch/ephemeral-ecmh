# Prover and KAT generator for the 14-byte binary variant (src/curvegen/select109.rs,
# src/curve/binary109/): per-namespace certificates plus hash-to-curve and sum
# vectors on each certified curve.
#
#   nix develop .#math -c sage sage/kat109.sage > tests/common/kats109.rs
#   nix fmt
#
# The field and group arithmetic is independent of the Rust code; the candidate
# derivation, the byte framing (src/hash.rs, src/curvegen/select109.rs) and the
# witness-selection rules match Rust by design. Rejection labels are the
# smallest prime factor of #E/2 without a bound (the Rust prover reproduces them
# given PARI's factors; src/curvegen/prove/gf2.rs); witnesses are the prover's
# hashed points (sage/kat_common.sage), so a run reproduces every byte.
load("sage/kat_common.sage")

MASK = 2^109 - 1

R.<t> = GF(2)[]
K.<z> = GF(2^109, modulus=t^109 + t^5 + t^4 + t^2 + 1)
TAG = b"ephemeral-ecmh/curve/gf2-109"


def encode(P):
    if P.is_zero():
        return 0
    x, y = P.xy()
    return gf2_to_int(x) | (ZZ(y.trace()) << 109)


def decode(E, x, sign):
    """The point with this x (Tr(x) = 1) and Tr(y) = sign, if any."""
    B = E.a6()
    if (B.sqrt() / x).trace() != 0:
        return None
    y = x * half_trace(x + 1 + B / x^2, 109)
    if ZZ(y.trace()) != sign:
        y += x
    return E(x, y)


def candidate(c):
    """x from bits 0..109, bit 0 fixed so that Tr(x) = 1; sign = bit 109."""
    x = gf2_from_int(K, c & MASK)
    if x.trace() == 0:
        x += 1
    return x, (c >> 109) & 1


def hash_(E, tag, salt, msg):
    ctr = 0
    while True:
        for c in halves(salted(tag, salt, msg, ctr)):
            P = decode(E, *candidate(c))
            if P is not None:
                return P
        ctr += 1


def prove(seed):
    rejections, j = [], 0
    while True:
        v = digest_half(TAG, seed, j) & MASK
        assert v != 0
        E = EllipticCurve(K, [1, 1, 0, 0, gf2_from_int(K, v)])
        n = E.order()
        r = n // 2
        if r.is_pseudoprime():
            assert mov_ok(2^109, r)
            # the verifier's own check: its hashed point has order r
            P = hash_(E, b"ephemeral-ecmh/cert-point", seed, int(j).to_bytes(4, "little"))
            assert not P.is_zero() and (r * P).is_zero()
            return E, j, r, rejections
        # hashed points have Tr(x) = 1 and lie in 2E, of odd order n/2
        l = small_factor(r)
        P = point_of_order(n // 2, l, lambda i: hash_(E, TAG_WITNESS, seed, witness_msg(j, i)))
        rejections.append((l, encode(P)))
        j += 1


# --- emit ------------------------------------------------------------------
out = sys.stdout
write_struct(
    out,
    "kat109.sage",
    "tests/common/kats109.rs",
    rejection_doc=["    /// (l, 14-byte encoding of P_j as a little-endian integer)"],
)
out.write("\npub const GF2_109_CERTS: &[CertKat] = &[\n")
for seed in SEEDS:
    E, i, r, rej = prove(seed)
    Ps = [hash_(E, b"kat", seed, bytes([k])) for k in range(NMSG)]
    write_cert(out, seed, i, r, rej, [encode(P) for P in Ps], encode(sum(Ps)), 28)
    print("seed %s.. index %d, %d rejections, r = 2^%.2f" % (seed[:4].hex(), i, len(rej), float(log(r, 2))),
          file=sys.stderr)
out.write("];\n")
