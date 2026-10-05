# Pornin's deterministic map (src/curve/binary/map.rs) from its definition,
# for each family's constants (src/curve/binary*/map.rs): the trace forms
# that make the constants work, and known answers on a fixed curve.
#
#   nix develop .#math -c sage sage/pornin_map.sage 127
#   nix develop .#math -c sage sage/pornin_map.sage 109
#   nix develop .#math -c sage sage/pornin_map.sage 122
#   nix develop .#math -c sage sage/pornin_map.sage 122-gls
#
# The report goes to stderr; stdout is the `VECTORS` table of the family's
# map tests, (input, x, y) of the affine output, as canonical integers.
#
# The point is decoded from w as Pornin's codec does (ePrint 2022/1325,
# §4.3; src/curve/binary/wcodec.rs), not through the extended coordinates
# the Rust map returns: x is the root of x^2 + d x + b with Tr(x) = 1, and
# y = x (x + lambda) with lambda = w^2 + 1 + a.
import sys

say = lambda *a: print(*a, file=sys.stderr)
which = sys.argv[1] if len(sys.argv) > 1 else "109"
assert which in ("127", "109", "122", "122-gls")
R.<t> = GF(2)[]


def roots(c):
    """The roots of X^2 + X = c."""
    return (X^2 + X + c).roots(multiplicities=False)


def pornin(a, B, c, sign, k, root_bit):
    """The map's affine output for the forced c and sign."""
    b = B.sqrt()
    ms = [c, c + k, c + c^2 / k]
    assert all(m.trace() == a.trace() == 1 for m in ms)
    es = [b / m for m in ms]
    assert sum(es) == 0
    i = min(j for j in range(3) if es[j].trace() == 0)
    d = ms[i].sqrt()
    (w,) = [r for r in roots(d + a) if root_bit(r) == sign]
    lam = w^2 + 1 + a
    (x,) = [r for r in (Y^2 + d * Y + b).roots(multiplicities=False) if r.trace() == 1]
    y = x * (x + lam)
    assert y^2 + x * y == x^3 + a * x^2 + B
    return x, y


# Inputs: tests/kat.rs's map_inputs, the first 8.
INPUTS = [((i + 1) * 0x9E3779B97F4A7C15F39CC0605CEDC835) % 2^128 for i in range(8)]
B_SEED = 0x5F3C8A712E94D0B61C47A389E5F26D1B

if which in ("127", "109"):
    m, f, sign_bit = {
        "127": (127, t^127 + t^63 + 1, 127),
        "109": (109, t^109 + t^5 + t^4 + t^2 + 1, 109),
    }[which]
    K.<z> = GF(2^m, modulus=f)
    P.<X> = K[]
    Q.<Y> = K[]
    to_int = lambda v: ZZ(v.polynomial().change_ring(ZZ)(2))
    from_int = lambda n: K(ZZ(n % 2^m).digits(2))

    # k = z^2 needs Tr(z^2) = Tr(z) = 0; then Tr(c) = Tr(a) = 1 and
    # Tr(c^2/z^2) = Tr(c/z) = 0 are each the parity of a few bits of c.
    k = z^2
    assert k.trace() == 0
    tr = [i for i in range(m) if (z^i).trace() == 1]
    trz = [i for i in range(m) if (z^i / z).trace() == 1]
    say("GF(2^%d): Tr(c) = sum of c_i, i in %s; Tr(c/z) = sum of c_i, i in %s" % (m, tr, trz))
    say("  z^-1 = %s, of trace %d" % (1 / z, (1 / z).trace()))
    if which == "127":
        assert tr == [0] and trz == [1]
    else:
        assert tr == [0, 105, 107] and trz == [1, 106, 108]

    def force(v):
        # bit 0 fixes Tr(c), bit 1 fixes Tr(c/z): each is in one form only
        v = v % 2^m
        bit = lambda i: (v >> i) & 1
        c0 = 1 ^^ sum(bit(i) for i in tr if i != 0) % 2
        c1 = sum(bit(i) for i in trz if i != 1) % 2
        return from_int((v >> 2 << 2) | c0 | c1 << 1)

    # Tr(1) = m mod 2 = 1: the trace tells w from w + 1
    root_bit = lambda w: ZZ(w.trace())
    a = K(1)
    B = from_int(B_SEED)
    say("  B = %#x" % to_int(B))
    rows = []
    for v in INPUTS:
        x, y = pornin(a, B, force(v), (v >> sign_bit) & 1, k, root_bit)
        rows.append((v, to_int(x), to_int(y)))

if which in ("122", "122-gls"):
    # GF(2^122) as one field, with the tower's z and u found in it: any
    # choice of the roots gives the same integers, as the map commutes
    # with the automorphisms that relate the choices.
    f = t^61 + t^23 + t^15 + t^5 + 1
    K = GF(2^122, "g")
    P.<X> = K[]
    Q.<Y> = K[]
    z = f.change_ring(K).roots(multiplicities=False)[0]
    u = (X^2 + X + 1).roots(multiplicities=False)[0]
    basis = [z^i for i in range(61)] + [z^i * u for i in range(61)]
    Mb = matrix(GF(2), [b._vector_() for b in basis])

    def to_int(v):
        bits = Mb.solve_left(v._vector_())
        return sum(ZZ(bits[i]) << i for i in range(61)) + sum(ZZ(bits[61 + i]) << (64 + i) for i in range(61))

    def from_int(n):
        sub = lambda w: sum(z^i for i in range(61) if (w >> i) & 1)
        return sub(n % 2^61) + sub((n >> 64) % 2^61) * u

    # Tr_122(a0 + a1 u) = Tr_61(a1), and Tr_61(z^i) = 1 for i < 61 only at
    # i = 0 (sage/gf2_122.sage): Tr(c) is bit 64 of the canonical c.
    tr = [i for i in range(128) if i % 64 < 61 and from_int(2^i).trace() == 1]
    # k = z^2 lies in GF(2^61), so Tr(k) = 0, and c/z divides both halves:
    # z^-1 = z^60 + z^22 + z^14 + z^4 has Tr_61 = 0, so Tr(c/z) is bit 65.
    trz = [i for i in range(128) if i % 64 < 61 and (from_int(2^i) / z).trace() == 1]
    k = z^2
    assert k.trace() == 0 and (1 / z).trace() == 0
    say("GF(2^122): Tr(c) = bits %s of c; Tr(c/z) = bits %s" % (tr, trz))
    assert tr == [64] and trz == [65]
    assert to_int(1 / z) == 2^60 + 2^22 + 2^14 + 2^4

    def force(v):
        # bit 64 = 1 and bit 65 = 0
        return from_int(((v >> 66) << 66) | (v % 2^64) | 2^64)

    # Tr(1) = 0 here, so the trace cannot tell w from w + 1: crrl takes
    # bit 0 of w0, which 1 flips.
    root_bit = lambda w: to_int(w) & 1
    a = u
    if which == "122":
        B = from_int(B_SEED)
    else:
        beta = from_int(B_SEED % 2^61)
        B = beta^4
    say("  B = %#x" % to_int(B))
    rows = []
    for v in INPUTS:
        x, y = pornin(a, B, force(v), (v >> 127) & 1, k, root_bit)
        rows.append((v, to_int(x), to_int(y)))

out = sys.stdout
out.write("    // sage/pornin_map.sage %s: (input, x, y) on B = %#x\n" % (which, to_int(B)))
out.write("    const VECTORS: [(u128, u128, u128); %d] = [\n" % len(rows))
for v, x, y in rows:
    out.write("        (%#x, %#x, %#x),\n" % (v, x, y))
out.write("    ];\n")
