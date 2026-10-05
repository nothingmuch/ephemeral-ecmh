# GF(2^122) = GF(2^61)[u]/(u^2 + u + 1) (src/field/gf2_122/): the base modulus,
# the tower's formulas, the base field's constant tables, and known answers.
#
#   nix develop .#math -c sage sage/gf2_122.sage tables > src/field/gf2_122/tables.rs
#   nix develop .#math -c sage sage/gf2_122.sage kats > src/field/gf2_122/kats.rs
#   nix fmt
#
# The report goes to stderr; stdout is the Rust file named by the argument.
import sys

say = lambda *a: print(*a, file=sys.stderr)
what = sys.argv[1] if len(sys.argv) > 1 else "tables"
assert what in ("tables", "kats")
R.<t> = GF(2)[]
m = 61

# --- base modulus ------------------------------------------------------------
# Swan: m = 61 = 5 mod 8, so every trinomial t^61 + t^k + 1 has an even number
# of irreducible factors.
assert m % 8 == 5
assert not any((t^m + t^k + 1).is_irreducible() for k in range(1, m))
say("irreducible trinomials of degree 61: none")
assert (t^m + t^5 + t^2 + t + 1).is_irreducible()
say("t^61 + t^5 + t^2 + t + 1: irreducible (the lowest pentanomial)")

# Elements live in one 64-bit word and are reduced modulo t^3 f, so a 128-bit
# product lo + hi t^64 folds onto rho = t^3 (f - t^61) with word-aligned
# shifts or carry-less multiplications: T = floor(hi rho / t^64), then
# lo + ((hi + T) rho mod t^64). That is exact while deg(T rho) < 64, i.e.
# 2 deg rho < 64 + 1, k3 <= 28, and then every pentanomial costs the same 4
# terms per fold. What differs is the trace and the square root:
#   - Tr(v) is the parity of v against {i < 64 : Tr(t^i) = 1};
#   - sqrt(v) = sqrt(v_even) + sqrt(t) sqrt(v_odd), one multiplication by
#     sqrt(t), which is a few shifts when sqrt(t) is sparse.
# With every middle exponent odd, t = t^62 + ... is a square times t^0 and
# sqrt(t) has 4 terms; Tr(t^i) = 0 for odd i < 61. Take the lowest such f.
def props(k3, k2, k1):
    K.<z> = GF(2^m, modulus=t^m + t^k3 + t^k2 + t^k1 + 1)
    tr = [i for i in range(64) if (z^i).trace() == 1]
    return tr, len(z.sqrt().polynomial().exponents())

pentas = [(k3, k2, k1) for k3 in range(3, 29) for k2 in range(2, k3) for k1 in range(1, k2)
          if (t^m + t^k3 + t^k2 + t^k1 + 1).is_irreducible()]
scored = sorted(pentas, key=lambda k: (len(props(*k)[0]) + props(*k)[1], k))
say("irreducible pentanomials with k3 <= 28: %d" % len(pentas))
for k in [(5, 2, 1)] + scored[:6]:
    tr, sq = props(*k)
    say("  t^61 + t^%d + t^%d + t^%d + 1: Tr bits %s, sqrt(t) has %d terms" % (k + (tr, sq)))
k3, k2, k1 = scored[0]
assert (k3, k2, k1) == (23, 15, 5)
f = t^m + t^k3 + t^k2 + t^k1 + 1
say("chosen: %s" % f)
F.<z> = GF(2^m, modulus=f)


def f_to_int(v):
    return ZZ(v.polynomial().change_ring(ZZ)(2))


def f_from_int(c):
    return F(ZZ(c).digits(2))


def half_trace(v):
    return sum(v^(2^(2 * i)) for i in range(31))


rho = 2^3 * (2^k3 + 2^k2 + 2^k1 + 1)
tr_bits = [i for i in range(64) if (z^i).trace() == 1]
assert tr_bits == [0, 61]
say("Tr(z^i) = 1 for i < 64 exactly at", tr_bits, "(Tr v = v_0 on canonical v)")
sqrt_z = z.sqrt()
assert sqrt_z == z^31 + z^12 + z^8 + z^3 and sqrt_z^2 == z
say("sqrt(z) =", sqrt_z)
assert half_trace(F(1)) == 1  # 31 terms

# --- the tower -----------------------------------------------------------------
# u^2 + u + 1 has no root in GF(2^61): a root would solve x^2 + x = 1, which
# needs Tr(1) = m mod 2 = 0.
P.<U> = F[]
assert (U^2 + U + 1).is_irreducible()
L.<u> = P.quotient(U^2 + U + 1)


def l_parts(a):
    c = a.lift().padded_list(2)
    return c[0], c[1]


def l_to_int(a):
    a0, a1 = l_parts(a)
    return f_to_int(a0) | f_to_int(a1) << 64


def frob(a, k):
    for _ in range(k):
        a = a^2
    return a


def l_trace(a):
    s = sum(frob(a, i) for i in range(2 * m))
    a0, a1 = l_parts(s)
    assert a1 == 0 and a0 in (0, 1)
    return ZZ(a0 == 1)


set_random_seed(122)
for _ in range(20):
    a0, a1, b0, b1 = [F.random_element() for _ in range(4)]
    a, b = a0 + a1 * u, b0 + b1 * u
    # Karatsuba: 3 base products
    c0, c1 = l_parts(a * b)
    assert c0 == a0 * b0 + a1 * b1 and c1 == (a0 + a1) * (b0 + b1) + a0 * b0
    # squaring, square root
    assert l_parts(a^2) == ((a0 + a1)^2, a1^2)
    assert l_parts(frob(a, 2 * m - 1)) == ((a0 + a1).sqrt(), a1.sqrt())
    # the conjugate is the 2^61-power: u -> u + 1; the norm is in GF(2^61)
    assert frob(a, m) == (a0 + a1) + a1 * u
    assert a * frob(a, m) == a0^2 + a0 * a1 + a1^2
    # Tr_122 = Tr_61(Tr_{122/61}), Tr_{122/61}(a) = a + a^(2^61) = a1
    assert a + frob(a, m) == a1
    assert l_trace(a) == ZZ(a1.trace())
    # z^2 + z = d: z1 = H(d1) solves z1^2 + z1 = d1 + Tr(d1); flipping z1 by 1
    # flips Tr(z1) (Tr_61(1) = 1), so pick it with Tr(z1) = Tr(d0), and then
    # z0 = H(d0 + z1^2) solves z0^2 + z0 = d0 + z1^2. So z^2 + z = d + u Tr(d)
    # always, and z^2 + z = d is solvable exactly when Tr_122(d) = Tr(d1) = 0.
    x1 = half_trace(a1)
    if x1.trace() != a0.trace():
        x1 += 1
    x0 = half_trace(a0 + x1^2)
    x = x0 + x1 * u
    assert x^2 + x == a + u * l_trace(a)
assert l_trace(u) == 1
say("tower: Karatsuba, square, sqrt, norm, trace and z^2 + z = d formulas agree with Sage")

out = sys.stdout


def tables():
    out.write("// @generated by sage/gf2_122.sage; do not edit. Regenerate from the repository root:\n")
    out.write("//   nix develop .#math -c sage sage/gf2_122.sage tables > src/field/gf2_122/tables.rs && nix fmt\n\n")
    out.write("/// rho = z^3 (z^%d + z^%d + z^%d + 1) = z^64 mod z^3 f\n" % (k3, k2, k1))
    out.write("pub const RHO: u64 = 0x%x;\n\n" % rho)
    out.write("/// H(z^i mod f) for i in 0..64: the half-trace of each word bit\n")
    out.write("pub const HALFTRACE: [u64; 64] = [\n")
    for i in range(64):
        out.write("    0x%016x,\n" % f_to_int(half_trace(z^i)))
    out.write("];\n")


# --- known answers -------------------------------------------------------------
def w(fmt, rows):
    for r in rows:
        out.write("    (" + ", ".join(fmt % v for v in r) + "),\n")


def kats():
    set_random_seed(1220)
    N = 8
    fs = [F.random_element() for _ in range(2 * N)] + [F(0), F(1), z, z^60]
    ls = [L(F.random_element() + F.random_element() * u) for _ in range(2 * N)] + [L(0), L(1), u, u + 1]
    out.write("// @generated by sage/gf2_122.sage kats; do not edit. Regenerate from the repository root:\n")
    out.write("//   nix develop .#math -c sage sage/gf2_122.sage kats > src/field/gf2_122/kats.rs && nix fmt\n")
    out.write("// GF(2^61) values are canonical; GF(2^122) ones are a0 | a1 << 64.\n\n")
    out.write("/// (a, b, a b) in GF(2^61)\npub const MUL61: &[(u64, u64, u64)] = &[\n")
    w("0x%x", [(f_to_int(a), f_to_int(b), f_to_int(a * b)) for a, b in zip(fs[::2], fs[1::2])])
    out.write("];\n\n/// (a, 1/a, sqrt(a), H(a), Tr(a)) in GF(2^61); 1/0 = 0\n")
    out.write("pub const UNARY61: &[(u64, u64, u64, u64, u32)] = &[\n")
    w("0x%x", [(f_to_int(a), f_to_int(a^-1 if a else a), f_to_int(a.sqrt()), f_to_int(half_trace(a)), ZZ(a.trace()))
               for a in fs])
    out.write("];\n\n/// (a, b, a b) in GF(2^122)\npub const MUL: &[(u128, u128, u128)] = &[\n")
    w("0x%x", [(l_to_int(a), l_to_int(b), l_to_int(a * b)) for a, b in zip(ls[::2], ls[1::2])])
    out.write("];\n\n/// (a, a^2, 1/a, sqrt(a), Tr(a)) in GF(2^122); 1/0 = 0\n")
    out.write("pub const UNARY: &[(u128, u128, u128, u128, u32)] = &[\n")
    w("0x%x", [(l_to_int(a), l_to_int(a^2), l_to_int(a^-1 if a else a), l_to_int(frob(a, 2 * m - 1)), l_trace(a))
               for a in ls])
    out.write("];\n")


tables() if what == "tables" else kats()
