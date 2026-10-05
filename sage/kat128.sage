# Independent certificates and known-answer vectors for E_d/<(0,-1)>, where
# E_d: -u^2 + v^2 = 1 + d*u^2*v^2 over p = 2^128 - 275.
#
# Run as Python through Sage (cypari2 supplies point counting and primality):
#   sage -python sage/kat128.sage --self-test
#   sage -python sage/kat128.sage > tests/common/kats128.rs
#   nix fmt
#
# The reference uses integer modular arithmetic and the affine Edwards law.
# It shares with Rust the specified byte framing, the candidate derivation and
# the witness-selection rules. PARI independently
# checks the birational Weierstrass model and scalar multiplication.
import hashlib
import math
import sys

from cypari2 import Pari


p = 2**128 - 275
HASSE = math.isqrt(4 * p)
R_LO = (p + 1 - HASSE + 3) // 4
R_HI = (p + 1 + HASSE) // 4
TAG = b"ephemeral-ecmh/curve/fp128"
TAG_CERT = b"ephemeral-ecmh/cert-point"
TAG_WITNESS = b"ephemeral-ecmh/prove/witness"
NMSG = 8
SEEDS = [
    bytes.fromhex("000000000019d6689c085ae165831e934ff763ae46a2a6c172b3f1b60a8ce26f"),
    bytes.fromhex("00000000839a8e6886ab5951d76f411475428afc90947ee320161bbf18eb6048"),
    hashlib.sha256(b"ephemeral-ecmh/kat/2").digest(),
    hashlib.sha256(b"ephemeral-ecmh/kat/3").digest(),
]
pari = Pari()
pari.allocatemem(256 * 1024 * 1024, silent=True)


def salted(tag, salt, msg, ctr):
    t = hashlib.sha256(tag).digest()
    return hashlib.sha256(t + t + salt + bytes(32) + msg + int(ctr).to_bytes(4, "little")).digest()


def halves(digest):
    return [int.from_bytes(digest[:16], "little"), int.from_bytes(digest[16:], "little")]


def square_root(a, q):
    """A square root for q = 5 mod 8; None for a nonsquare."""
    assert q % 8 == 5
    a %= q
    root = pow(a, (q + 3) // 8, q)
    if root * root % q != a:
        root = root * pow(2, (q - 1) // 4, q) % q
    return root if root * root % q == a else None


def on_curve(P, d, q=p):
    u, v = P
    return (-u * u + v * v - 1 - d * u * u * v * v) % q == 0


def add(P, Q, d, q=p):
    u, v = P
    x, y = Q
    t = d * u * v * x * y % q
    return ((u * y + v * x) * pow((1 + t) % q, -1, q) % q,
            (v * y + u * x) * pow((1 - t) % q, -1, q) % q)


def mul(n, P, d, q=p):
    R = (0, 1)
    while n:
        if n & 1:
            R = add(R, P, d, q)
        P = add(P, P, d, q)
        n >>= 1
    return R


def encode(P, q=p):
    u, v = P
    return -u % q if (v & 1) or (v == 0 and u & 1) else u


def decode(c, d, q=p):
    if not 0 <= c < q:
        return None
    v2 = (1 + c * c) * pow((1 - d * c * c) % q, -1, q) % q
    v = square_root(v2, q)
    if v is None or (v == 0 and c & 1):
        return None
    if v & 1:
        v = -v % q
    return c, v


def hash_point(d, seed, msg, tag=b"kat"):
    for ctr in range(2**32):
        for c in halves(salted(tag, seed, msg, ctr)):
            P = decode(c, d)
            if P is not None:
                return P
    raise RuntimeError("hash counter exhausted")


def weierstrass(d, q=p):
    """Montgomery U=(1+v)/(1-v), V=U/u; X=U/B, Y=V/B.

    B=-4/(1+d), A=2*(1-d)/(1+d). Thus the Weierstrass equation is
    Y^2 = X^3 + (d-1)/2*X^2 + (d+1)^2/16*X.
    """
    a2 = (d - 1) * pow(2, -1, q) % q
    a4 = (d + 1)**2 * pow(16, -1, q) % q
    return pari.ellinit([0, a2, 0, a4, 0], q)


def to_weierstrass(P, d, q=p):
    u, v = P
    if P == (0, 1):
        return pari([0])
    if u == 0:
        assert v == q - 1
        return pari([0, 0])
    x = -(1 + d) * (1 + v) * pow(4 * (1 - v) % q, -1, q) % q
    return pari([x, x * pow(u, -1, q) % q])


def small_factor(n):
    for ell in pari.primes([3, 2**16]):
        ell = int(ell)
        if n % ell == 0:
            return ell
    return int(pari.factor(n)[0][0])


def witness(d, n, ell, seed, index):
    """A deterministic raw point of order ell, for ell prime or ell=8."""
    base = 2 if ell == 8 else ell
    m = n
    while m % base == 0:
        m //= base
    for attempt in range(2**32):
        msg = index.to_bytes(4, "little") + attempt.to_bytes(4, "little")
        P = mul(m, hash_point(d, seed, msg, TAG_WITNESS), d)
        if (mul(4, P, d) if ell == 8 else P) == (0, 1):
            continue
        while mul(ell, P, d) != (0, 1):
            P = mul(base, P, d)
        c = encode(P)
        Q = decode(c, d)
        assert on_curve(Q, d) and encode(Q) == c
        if ell == 8:
            assert mul(2, Q, d)[0] != 0 and mul(4, Q, d)[0] == 0
        else:
            assert 3 <= ell < R_LO and Q[0] != 0 and mul(ell, Q, d)[0] == 0
        E = weierstrass(d)
        assert int(pari.ellorder(E, to_weierstrass(P, d), n)) == ell
        return ell, c
    raise RuntimeError("witness counter exhausted")


def prove(seed):
    rejections = []
    for j in range(2**32):
        d = halves(salted(TAG, seed, b"", j))[0] % p
        if pow(d, (p - 1) // 2, p) != p - 1:
            continue
        E = weierstrass(d)
        n = int(pari.ellcard(E))
        assert n % 4 == 0 and abs(n - p - 1) <= HASSE
        r = n // 4
        if n % 8:
            if pari.isprime(r):
                assert R_LO <= r <= R_HI
                assert pari.znorder(pari.Mod(p, r)) > 2**20  # embedding degree, as EMBEDDING_MIN
                Q = mul(2, hash_point(d, seed, j.to_bytes(4, "little"), TAG_CERT), d)
                assert Q[0] != 0 and mul(r, Q, d)[0] == 0
                assert pari.ellmul(E, to_weierstrass(Q, d), r) in (pari([0]), pari([0, 0]))
                return d, j, r, rejections
            ell = small_factor(r)
        else:
            ell = 8
        rejections.append(witness(d, n, ell, seed, j))
        print("candidate %d rejected by %d" % (j, ell), file=sys.stderr, flush=True)
    raise RuntimeError("candidate counter exhausted")


def self_test():
    """Exhaustively compare the codec with enumerated two-element cosets."""
    curves = pairs = 0
    for q in (5, 13, 29, 37):
        for d in range(1, q):
            if pow(d, (q - 1) // 2, q) != q - 1:
                continue
            Ps = [(u, v) for u in range(q) for v in range(q) if on_curve((u, v), d, q)]
            E = weierstrass(d, q)
            assert int(pari.ellcard(E)) == len(Ps)
            cosets = {P: frozenset((P, (-P[0] % q, -P[1] % q))) for P in Ps}
            codes = {}
            for P in Ps:
                c = encode(P, q)
                assert decode(c, d, q) in cosets[P]
                assert (c == 0) == (P[0] == 0)
                assert pari.ellisoncurve(E, to_weierstrass(P, d, q))
                assert codes.setdefault(c, cosets[P]) == cosets[P]
            assert len(codes) * 2 == len(Ps)
            for c in range(-1, q + 1):
                assert (decode(c, d, q) is not None) == (c in codes)
            for P in Ps:
                for Q in Ps:
                    A, B = decode(encode(P, q), d, q), decode(encode(Q, q), d, q)
                    R = add(P, Q, d, q)
                    assert on_curve(R, d, q)
                    assert to_weierstrass(R, d, q) == pari.elladd(E, to_weierstrass(P, d, q), to_weierstrass(Q, d, q))
                    assert encode(R, q) == encode(add(A, B, d, q), q)
                    assert encode(add(P, (-Q[0] % q, Q[1]), d, q), q) == encode(add(A, (-B[0] % q, B[1]), d, q), q)
                    assert ((P[0] * Q[1] - Q[0] * P[1]) % q == 0) == (cosets[P] == cosets[Q])
                    pairs += 1
            curves += 1
    print("PASS: %d curves, %d ordered pairs; codec, signed closure, PARI group law" % (curves, pairs), file=sys.stderr)


def emit():
    assert pari.isprime(p) and HASSE == 2**65 - 1
    assert R_LO == 2**126 - 2**63 - 68 and R_HI == 2**126 + 2**63 - 69
    print("// @generated by sage/kat128.sage; do not edit. Regenerate from the repository root:")
    print("//   nix develop .#math -c sage -python sage/kat128.sage > tests/common/kats128.rs && nix fmt")
    print("#![allow(dead_code)]\n")
    print("/// As `kats::CertKat`, with 16-byte two-torsion quotient encodings.")
    print("pub struct CertKat {\n    pub seed: [u8; 32],\n    pub index: u32,\n    pub r: u128,")
    print("    pub rejections: &'static [(u128, u128)],")
    print('    /// hash_to_curve(Salted::new(b"kat", seed), [i]) for i in 0..8, encoded')
    print("    pub hashes: [u128; 8],\n    /// Encoded sum of `hashes`.\n    pub sum: u128,\n}\n")
    print("pub const FP128_CERTS: &[CertKat] = &[")
    overflow = False
    for seed in SEEDS:
        d, index, r, rej = prove(seed)
        Ps = [hash_point(d, seed, bytes([k])) for k in range(NMSG)]
        total = (0, 1)
        E = weierstrass(d)
        pari_total = pari([0])
        for P in Ps:
            assert on_curve(P, d) and decode(encode(P), d) == P
            total = add(total, P, d)
            pari_total = pari.elladd(E, pari_total, to_weierstrass(P, d))
        assert to_weierstrass(total, d) == pari_total
        overflow |= 4 * r >= 2**128
        print("    CertKat {\n        seed: %s," % list(seed))
        print("        index: %d,\n        r: 0x%x,\n        rejections: &[" % (index, r))
        for ell, c in rej:
            print("            (0x%x, 0x%032x)," % (ell, c))
        print("        ],\n        hashes: [")
        for P in Ps:
            print("            0x%032x," % encode(P))
        print("        ],\n        sum: 0x%032x,\n    }," % encode(total))
        print("FP128: seed %s.. index %d, %d rejections, #E %s 2^128" % (seed[:4].hex(), index, len(rej), ">=" if 4 * r >= 2**128 else "<"), file=sys.stderr, flush=True)
    print("];\n")
    assert overflow, "include an additional seed to exercise #E >= 2^128"


if __name__ == "__main__":
    if sys.argv[1:] == ["--self-test"]:
        self_test()
    elif not sys.argv[1:]:
        emit()
    else:
        raise SystemExit("usage: sage -python sage/kat128.sage [--self-test]")
