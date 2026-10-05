# Prover and KAT generator for the Rust implementation: per-namespace curve
# certificates (src/curvegen/select.rs) plus hash-to-curve and sum vectors on
# each certified curve.
#
#   nix develop .#math -c sage sage/kat.sage > tests/common/kats.rs
#   nix fmt
#
# The field and group arithmetic is Sage's, independent of the Rust code; the
# candidate derivation, the byte framing (src/hash.rs, src/curvegen/select.rs)
# and the witness-selection rules match Rust by design. Rejection witnesses
# are the Rust prover's hashed points (sage/kat_common.sage), so a run
# reproduces every byte.
load("sage/kat_common.sage")

MASK = 2^127 - 1

# --- GF(2^127) -------------------------------------------------------------
R.<t> = GF(2)[]
K.<z> = GF(2^127, modulus=t^127 + t^63 + 1)
TAG_GF2 = b"ephemeral-ecmh/curve/gf2"


def gf2_encode(P):
    if P.is_zero():
        return 0
    x, y = P.xy()
    return gf2_to_int(x) | (ZZ(y.polynomial().padded_list(127)[0]) << 127)


def gf2_decode(E, c):
    B = E.a6()
    b = B.sqrt()
    if c == 0:
        return E(0)
    if c & 1 == 0:
        return None
    x, sign = gf2_from_int(K, c & MASK), c >> 127
    if (b / x).trace() != 0:
        return None
    y = x * half_trace(x + 1 + B / x^2, 127)
    if ZZ(y.polynomial().padded_list(127)[0]) != sign:
        y += x
    return E(x, y)


def gf2_hash(E, salt, msg, tag=b"kat"):
    ctr = 0
    while True:
        for c in halves(salted(tag, salt, msg, ctr)):
            P = gf2_decode(E, c | 1)
            if P is not None:
                return P
        ctr += 1


def gf2_prove(seed):
    rejections, j = [], 0
    while True:
        v = digest_half(TAG_GF2, seed, j) & MASK
        assert v != 0
        E = EllipticCurve(K, [1, 1, 0, 0, gf2_from_int(K, v)])
        n = E.order()
        r = n // 2
        if r.is_pseudoprime():
            assert mov_ok(2^127, r)
            return E, j, r, rejections
        # hashed points have Tr(x) = 1 and lie in 2E, of odd order n/2
        l = small_factor(r)
        P = point_of_order(n // 2, l, lambda i: gf2_hash(E, seed, witness_msg(j, i), TAG_WITNESS))
        rejections.append((l, gf2_encode(P)))
        j += 1


# --- emit ------------------------------------------------------------------
out = sys.stdout
write_struct(out, "kat.sage", "tests/common/kats.rs")
out.write("\n")


def emit(name, prove, hash_, encode):
    out.write("pub const %s: &[CertKat] = &[\n" % name)
    for seed in SEEDS:
        E, i, r, rej = prove(seed)
        Ps = [hash_(E, seed, bytes([k])) for k in range(NMSG)]
        write_cert(out, seed, i, r, rej, [encode(P) for P in Ps], encode(sum(Ps)), 32)
        print("%s: seed %s.. index %d, %d rejections" % (name, seed[:4].hex(), i, len(rej)), file=sys.stderr)
    out.write("];\n\n")


emit("GF2_127_CERTS", gf2_prove, gf2_hash, gf2_encode)
