# Prover and KAT generator for the GF(2^122) families (src/curvegen/select122.rs,
# src/curve/binary122/): per-namespace certificates, dense and GLS, plus
# hash-to-curve and sum vectors on each certified curve.
#
#   nix develop .#math -c sage sage/kat122.sage > tests/common/kats122.rs
#   nix fmt
#
# The field and group arithmetic is independent of the Rust code; the candidate
# derivation, the byte framing (src/hash.rs, src/curvegen/select122.rs), the
# tower's integer form (a0 | a1 << 64 for a0 + a1 u, which src/field/gf2_122
# fixes) and the witness-selection rules match Rust by design. Arithmetic is
# Sage's GF(2^122), with GF(2^61) and u embedded; the certified orders come from
# E.order() and, for GLS, are checked against (q - 1)^2 + t^2 over the subfield.
# Each certified curve's magic numbers (sage/ghs122.sage) go to stderr.
#
# Rejection labels are the smallest prime factor of #E/2 below 2^16, else #E
# itself (the even-order marker); witnesses are the prover's hashed points
# (sage/kat_common.sage), so a run reproduces every byte. A j that is no
# candidate takes a (0, 0) entry.
load("sage/kat_common.sage")

say = lambda *a: print(*a, file=sys.stderr)
MASK61 = 2^61 - 1
MASK122 = MASK61 | MASK61 << 64
Q = 2^122
q = 2^61

R.<t> = GF(2)[]
F61.<w> = GF(2^61, modulus=t^61 + t^23 + t^15 + t^5 + 1)
K.<g> = GF(2^122)
# the tower inside K: a root of the base modulus, and u with u^2 + u + 1 = 0
W = sorted((t^61 + t^23 + t^15 + t^5 + 1).change_ring(K).roots(multiplicities=False))[0]
U = sorted((t^2 + t + 1).change_ring(K).roots(multiplicities=False))[0]
BASIS = [W^i for i in range(61)] + [W^i * U for i in range(61)]
V = K.vector_space(map=True)[2]
TO_TOWER = matrix(GF(2), [V(b) for b in BASIS]).inverse()


def f61(c):
    return F61(ZZ(c).digits(2))


def embed61(v):
    """GF(2^61) -> K along the tower's embedding."""
    return sum(ZZ(c) * W^i for i, c in enumerate(v.polynomial().padded_list(61)))


def k_from_int(c):
    """a0 | a1 << 64 -> a0 + a1 u."""
    a0, a1 = c & MASK61, (c >> 64) & MASK61
    return embed61(f61(a0)) + embed61(f61(a1)) * U


def k_to_int(v):
    bits = V(v) * TO_TOWER
    a0 = sum(ZZ(bits[i]) << i for i in range(61))
    a1 = sum(ZZ(bits[61 + i]) << i for i in range(61))
    return a0 | a1 << 64


assert k_to_int(U) == 1 << 64 and k_to_int(W) == 2 and k_to_int(K(1)) == 1
assert all(k_to_int(k_from_int(c)) == c for c in [ZZ.random_element(2^128) & MASK122 for _ in range(20)])


def encode(P):
    if P.is_zero():
        return 0
    x, y = P.xy()
    return k_to_int(x) | (ZZ(y.trace()) << 127)


def decode(E, x, sign):
    """The point with this x (Tr(x) = 1) and Tr(y) = sign, if any."""
    ys = [y for y in E.lift_x(x, all=True)]
    if not ys:
        return None
    (P,) = [P for P in ys if ZZ(P.xy()[1].trace()) == sign]
    return P


def candidate(c):
    """x from the bits of MASK122 with bit 64 set, so Tr(x) = 1; sign = bit 127."""
    x = k_from_int((c | 1 << 64) & MASK122)
    assert x.trace() == 1
    return x, c >> 127


def hash_(E, tag, salt, msg):
    ctr = 0
    while True:
        for c in halves(salted(tag, salt, msg, ctr)):
            P = decode(E, *candidate(c))
            if P is not None:
                return P
        ctr += 1


def curve(B):
    return EllipticCurve(K, [1, U, 0, 0, B])


# the two families' search orders: candidate j -> (B, beta in GF(2^61) or
# None), or None for a j that is no candidate
def dense_b(seed, j):
    """B = the first digest half, masked to the tower's bits; B in GF(4) is
    no candidate (a curve over GF(4), sage/ghs122.sage)."""
    v = digest_half(b"ephemeral-ecmh/curve/gf2-122", seed, j) & MASK122
    return None if v & (MASK122 - 1 - 2^64) == 0 else (k_from_int(v), None)


def gls_b(seed, j):
    """beta = the first digest half's low 61 bits, B = beta^4; beta in GF(2)
    is no candidate."""
    v = digest_half(b"ephemeral-ecmh/curve/gf2-122-gls", seed, j) & MASK61
    return None if v <= 1 else (embed61(f61(v))^4, f61(v))


def magic(B, l):
    """GHS's m for the descent to GF(2^l): dim Span{(1, sigma^i(sqrt B))}."""
    s, rows = B.sqrt(), []
    for _ in range(122 // l):
        rows.append([1] + list(V(s)))
        s = s^(2^l)
    return matrix(GF(2), rows).rank()


def prove(seed, family):
    rejections, j = [], 0
    while True:
        cand = (dense_b if family == "dense" else gls_b)(seed, j)
        if cand is None:
            rejections.append((0, 0))  # ignored by the verifier
            j += 1
            continue
        B, beta = cand
        E = curve(B)
        n = E.order()
        assert n % 4 == 2  # Tr(u) = 1: one point of order 2
        if beta is not None:
            # the quadratic twist of y^2 + xy = x^3 + beta^4 over GF(2^61)
            t = q + 1 - EllipticCurve(F61, [1, 0, 0, 0, beta^4]).order()
            assert n == (q - 1)^2 + t^2
        r = n // 2
        if r.is_pseudoprime():
            assert r.is_prime()
            k = Mod(Q, r).multiplicative_order()
            assert k > EMBEDDING_MIN
            P = hash_(E, b"ephemeral-ecmh/cert-point", seed, int(j).to_bytes(4, "little"))
            assert not P.is_zero() and (r * P).is_zero()
            return E, j, r, rejections, k
        # hashed points have Tr(x) = 1 and lie in 2E, of odd order n/2
        l = small_factor(r, full=False)
        if l is not None:
            P = point_of_order(n // 2, l, lambda i: hash_(E, TAG_WITNESS, seed, witness_msg(j, i)))
            rejections.append((l, encode(P)))
        else:
            # rejected by the order itself: n in the Hasse interval, n/2 not
            # prime; the prover's order witness is its first hashed point
            P = hash_(E, TAG_WITNESS, seed, witness_msg(j, 0))
            assert not P.is_zero() and (n * P).is_zero()
            rejections.append((n, encode(P)))
        j += 1


# --- emit ------------------------------------------------------------------
out = sys.stdout
write_struct(
    out,
    "kat122.sage",
    "tests/common/kats122.rs",
    rejection_doc=[
        "    /// (l, 16-byte encoding of P_j as a little-endian integer); (0, 0) for",
        "    /// a j that is no candidate",
    ],
    extra=[
        "    /// encoded 2 hashes[0] and hashes[0] - hashes[1]",
        "    pub double: u128,",
        "    pub diff: u128,",
    ],
)
for family, name in [("dense", "GF2_122_CERTS"), ("gls", "GF2_122_GLS_CERTS")]:
    out.write("\npub const %s: &[CertKat] = &[\n" % name)
    for seed in SEEDS:
        E, i, r, rej, k = prove(seed, family)
        Ps = [hash_(E, b"kat", seed, bytes([m])) for m in range(NMSG)]
        write_cert(
            out, seed, i, r, rej, [encode(P) for P in Ps], encode(sum(Ps)), 32,
            extra=[("double", encode(2 * Ps[0])), ("diff", encode(Ps[0] - Ps[1]))],
        )
        B = E.a6()
        say("%s seed %s.. index %d, %d rejections (%d by order), r = 2^%.2f, "
            "embedding degree 2^%.1f, magic numbers to GF(2), GF(4), GF(2^61): %d, %d, %d"
            % (family, seed[:4].hex(), i, len(rej), sum(1 for l, _ in rej if l % 2 == 0 and l > 0),
               float(log(r, 2)), float(log(k, 2)), magic(B, 1), magic(B, 2), magic(B, 61)))
    out.write("];\n")
