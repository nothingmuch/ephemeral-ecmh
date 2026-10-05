# What the Sage KAT generators share: the byte framing of src/hash.rs, the
# certificate search helpers, the seeds and the `CertKat` emitter. Loaded by
# sage/kat.sage, kat109.sage and kat122.sage, which are run from
# the repository root:
#
#   nix develop .#math -c sage sage/kat<family>.sage > tests/common/kats<family>.rs
#
# Each generator keeps what is its own: the field and its integer form, the
# curve model, the point codec, the candidate filter and the label rule.
# Rejection witnesses are the Rust prover's (src/curvegen/prove.rs): points
# hashed under TAG_WITNESS with message j || i, taken to the label's order
# by `point_of_order`. Nothing here draws random numbers, so a run
# reproduces every byte.
import hashlib, sys


# --- framing (src/hash.rs) ---------------------------------------------------
def salted(tag, salt, msg, ctr):
    t = hashlib.sha256(tag).digest()
    return hashlib.sha256(t + t + salt + bytes(32) + msg + int(ctr).to_bytes(4, "little")).digest()


def halves(d):
    return [int.from_bytes(d[:16], "little"), int.from_bytes(d[16:], "little")]


def digest_half(tag, seed, j):
    """The first half of the tagged digest of j: the raw material of candidate j."""
    return halves(salted(tag, seed, b"", j))[0]


# --- the search ----------------------------------------------------------------
def small_factor(n, full=True):
    """Smallest prime factor of odd composite n: by trial division below 2^16,
    then by factoring it, or None when `full` is off."""
    for l in primes(3, 2^16):
        if n % l == 0:
            return l
    return factor(n)[0][0] if full else None


EMBEDDING_MIN = 2^20


def mov_ok(q, r):
    """The embedding degree of r over F_q exceeds EMBEDDING_MIN (Rust's constant)."""
    return Mod(q, r).multiplicative_order() > EMBEDDING_MIN


TAG_WITNESS = b"ephemeral-ecmh/prove/witness"


def witness_msg(j, i):
    """The prover's message for the i-th hashed try on candidate j: j || i, u32 LE."""
    return int(j).to_bytes(4, "little") + int(i).to_bytes(4, "little")


def point_of_order(m, l, hash_):
    """A point of order exactly l (a prime, or 8) from hashed points of a
    group of order m, as the prover's `prime_power_point`: the i-th hashed
    point times m with its l-part removed has order a power of the prime
    under l; multiplying by that prime until l kills it leaves order l,
    unless the point was O or fell short of l, when the next i is tried.
    (m/l)*R alone is not enough: if E[l] = Z/l x Z/l, it is always O."""
    base = 2 if l == 8 else l
    k = m
    while k % base == 0:
        k //= base
    i = 0
    while True:
        P = k * hash_(i)
        i += 1
        if P.is_zero():
            continue
        while not (l * P).is_zero():
            P = base * P
        if not ((l // base) * P).is_zero():
            return P


def sqrt_canonical(a):
    """The square root Rust's prime fields take: a^((p + 1)/4), p = 3 mod 4; None for a nonsquare."""
    p = a.parent().characteristic()
    assert p % 4 == 3
    r = a^((p + 1) // 4)
    return r if r^2 == a else None


def edwards_order8(codec, E):
    """The prover's label-8 witness on an a = 1 Edwards curve, as
    `sieve::edwards_order8` finds it without a point count. On the
    Montgomery model y^2 = x (x^2 + a2 x + e^2), e = (1 - d)/4 = (1 - a2)/2,
    the points of order 8 have x = r (r + 1 + t) with r^2 = e and
    t^2 = 2r + 1: r the canonical root of e or its negative, the first
    with a rational t, and the point of sign 0 over that x. None when
    no such x exists, which is when 8 does not divide #E."""
    r = sqrt_canonical((1 - E.a2()) / 2)
    if r is None:
        return None
    for r in (r, -r):
        t = sqrt_canonical(2 * r + 1)
        if t is None:
            continue
        P = codec.decode(E, ZZ(r * (r + 1 + t)))
        if P is not None:
            return P
    return None


def gf2_from_int(K, c):
    """The element of K = GF(2^n) in polynomial basis whose coefficients are the bits of c."""
    return K(ZZ(c).digits(2))


def gf2_to_int(v):
    return ZZ(v.polynomial().change_ring(ZZ)(2))


def half_trace(a, n):
    """The half-trace of a in GF(2^n), n odd: a solution z of z^2 + z = a when
    Tr(a) = 0 (the other is z + 1); callers check the trace first."""
    return sum(a^(2^(2 * i)) for i in range((n + 1) // 2))


# --- emit ----------------------------------------------------------------------
SEEDS = [
    # block 0 and block 1 hashes, display (big-endian) byte order
    bytes.fromhex("000000000019d6689c085ae165831e934ff763ae46a2a6c172b3f1b60a8ce26f"),
    bytes.fromhex("00000000839a8e6886ab5951d76f411475428afc90947ee320161bbf18eb6048"),
    hashlib.sha256(b"ephemeral-ecmh/kat/2").digest(),
    hashlib.sha256(b"ephemeral-ecmh/kat/3").digest(),
]
NMSG = 8


def write_struct(out, script, output, struct_doc=(), rejection_doc=(), extra=()):
    """The `CertKat` header, `output` being the file the script writes.
    `struct_doc` lines precede the struct, `rejection_doc` lines the rejections
    field, and `extra` lines (doc and fields) follow `sum`. Lines are written
    as given, with the `/// ` and indentation."""
    out.write("// @generated by sage/%s; do not edit. Regenerate from the repository root:\n" % script)
    out.write("//   nix develop .#math -c sage sage/%s > %s && nix fmt\n" % (script, output))
    out.write("#![allow(dead_code)]\n\n")
    for line in struct_doc:
        out.write(line + "\n")
    out.write("pub struct CertKat {\n    pub seed: [u8; 32],\n    pub index: u32,\n    pub r: u128,\n")
    for line in rejection_doc:
        out.write(line + "\n")
    out.write("    pub rejections: &'static [(u128, u128)],\n")
    out.write("    /// hash_to_curve(Salted::new(b\"kat\", seed), [i]) for i in 0..%d, encoded\n" % NMSG)
    out.write("    pub hashes: [u128; %d],\n" % NMSG)
    out.write("    /// encoded sum of `hashes`\n    pub sum: u128,\n")
    for line in extra:
        out.write(line + "\n")
    out.write("}\n")


def write_cert(out, seed, i, r, rej, hashes, total, width, extra=()):
    """One `CertKat` entry; `hashes`, `total` and `extra` values are encoded
    integers, written as 0x with `width` hex digits."""
    width = int(width)
    out.write("    CertKat {\n        seed: %s,\n" % list(seed))
    out.write("        index: %d,\n        r: 0x%x,\n        rejections: &[\n" % (i, r))
    for l, c in rej:
        out.write("            (0x%x, 0x%0*x),\n" % (l, width, c))
    out.write("        ],\n        hashes: [\n")
    for h in hashes:
        out.write("            0x%0*x,\n" % (width, h))
    out.write("        ],\n        sum: 0x%0*x,\n" % (width, total))
    for name, value in extra:
        out.write("        %s: 0x%0*x,\n" % (name, width, value))
    out.write("    },\n")


# --- the odd prime fields --------------------------------------------------------
class PrimeCodec:
    """Points of a curve over GF(p), p = 2^bits - 1, as integers: x | sign << bits,
    O as p, the sign being y's parity (src/curve/encoding.rs). Masking to
    `bits` yields a value in 0..=p, so no range check beyond O's is needed."""

    def __init__(self, F, bits):
        assert F.characteristic() == 2^bits - 1
        self.F, self.p, self.bits, self.mask = F, F.characteristic(), bits, 2^bits - 1

    def encode(self, P):
        if P.is_zero():
            return self.p
        x, y = P.xy()
        return ZZ(x) | ((ZZ(y) & 1) << self.bits)

    def decode(self, E, c):
        v, sign = c & self.mask, c >> self.bits
        if v == self.p:
            return E(0) if sign == 0 else None
        x = self.F(v)
        rhs = x^3 + E.a2() * x^2 + E.a4() * x + E.a6()
        if not rhs.is_square():
            return None
        y = rhs.sqrt()
        if y == 0 and sign:
            return None
        if ZZ(y) & 1 != sign:
            y = -y
        return E(x, y)

    def hash(self, E, salt, msg, tag=b"kat"):
        """Each digest half, masked to bits + 1, decoded if it is not O's encoding."""
        ctr = 0
        while True:
            for c in halves(salted(tag, salt, msg, ctr)):
                c &= 2^(self.bits + 1) - 1
                if c & self.mask == self.p:
                    continue
                P = self.decode(E, c)
                if P is not None:
                    return P
            ctr += 1


def edwards_prove(codec, tag, seed):
    """Edwards d non-square -> y^2 = x^3 + a2 x^2 + a4 x; #E = 4r, r a pseudoprime.
    A candidate with 8 | #E is rejected under label 8 by the sieve's point
    (`edwards_order8`), else by the smallest prime factor of #E/4."""
    F, p, rejections, j = codec.F, codec.p, [], 0
    while True:
        d = F(digest_half(tag, seed, j) & codec.mask)
        if not d.is_square():
            E = EllipticCurve(F, [0, (1 + d) / 2, 0, (1 - d)^2 / 16, 0])
            n = E.order()
            P8 = edwards_order8(codec, E)
            assert (P8 is None) == (n % 8 != 0), "the order-8 rule disagrees with #E"
            if n % 8 == 0:
                assert (8 * P8).is_zero() and not (4 * P8).is_zero()
                rejections.append((8, codec.encode(P8)))
            else:
                r = n // 4
                if r.is_pseudoprime():
                    assert mov_ok(p, r)
                    return E, j, r, rejections
                l = small_factor(r)
                P = point_of_order(n, l, lambda i: codec.hash(E, seed, witness_msg(j, i), TAG_WITNESS))
                rejections.append((l, codec.encode(P)))
        j += 1


def weier_prove(codec, tag, seed):
    """y^2 = x^3 - 3x + b of pseudoprime order, not p. A candidate of even
    order is rejected under label 2, else by the smallest prime factor of #E."""
    F, p, rejections, j = codec.F, codec.p, [], 0
    while True:
        b = F(digest_half(tag, seed, j) & codec.mask)
        if b not in (0, 2, -2):
            E = EllipticCurve(F, [-3, b])
            n = E.order()
            if n.is_pseudoprime() and n != p:
                assert mov_ok(p, n)
                return E, j, n, rejections
            l = 2 if n % 2 == 0 else small_factor(n)
            P = point_of_order(n, l, lambda i: codec.hash(E, seed, witness_msg(j, i), TAG_WITNESS))
            rejections.append((l, codec.encode(P)))
        j += 1
