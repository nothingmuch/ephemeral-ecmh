# GHS / Weil descent check for random curves over GF(2^m), m prime.
#
#   nix develop .#math -c sage sage/ghs.sage
#
# GHS check for curves y^2 + xy = x^3 + x^2 + B over GF(2^m), m prime.
# Only subfield is GF(2), so the descent is to GF(2) with n = m, genus 2^(mb-1)
# (or 2^(mb-1)-1), mb = magic number = dim span_F2{(1, sqrt(B)^(2^i))}.
# sigma-invariant subspaces have dims built from the factors of x^m - 1 over F2:
# 1 (for x - 1) and d = ord_m(2) (for each factor of Phi_m), so mb is 1, d, d+1, ...
R.<x> = GF(2)[]
print("m    ord_m(2)=d  #deg-d factors  min mb for B not in GF(2)  -> min genus")
for m in [89, 97, 101, 103, 107, 109, 113, 127]:
    d = Mod(2, m).multiplicative_order()
    print("%-4d %-11d %-15d %-26s 2^%d" % (m, d, (m - 1) // d, "%d or %d" % (d, d + 1), d - 1))

def magic(K, B):
    s = B.sqrt(); rows = []; c = s
    for _ in range(K.degree()):
        rows.append([1] + c.polynomial().padded_list(K.degree())); c = c^2
    return matrix(GF(2), rows).rank()

m = 127
K.<z> = GF(2^m, modulus=x^127 + x^63 + 1)
V, from_V, to_V = K.vector_space(map=True)
S = matrix(GF(2), [to_V(from_V(e)^2) for e in V.basis()]).transpose()  # Frobenius
facs = [f for f, _ in (x^m - 1).factor() if f.degree() > 1]
print("\nm = 127: x^127 - 1 = (x+1) * %d irreducible factors of degree %s" % (len(facs), sorted(set(f.degree() for f in facs))))

# weak classes: sqrt(B) in V_i or 1 + V_i (V_i = ker f_i(sigma), dim 7); one Frobenius orbit each
from collections import Counter
weak = []
for f in facs:
    v = from_V(f(S).right_kernel().basis()[0])
    for s in (v, 1 + v):
        B = s^2
        n = EllipticCurve(K, [1, 1, 0, 0, B]).order()
        weak.append((magic(K, B), n))
print("weak classes: %d, magic numbers %s" % (len(weak), dict(Counter(mb for mb, _ in weak))))
print("  of which #E = 2*prime: %d" % sum(1 for _, n in weak if n % 4 == 2 and (n // 2).is_pseudoprime()))
W = set(n for _, n in weak)
print("  distinct orders (isogeny classes to exclude): %d" % len(W))

# the candidate beta = z^43 and random B
n43 = EllipticCurve(K, [1, 1, 0, 0, (z^43)^4]).order()
print("\nbeta = z^43: magic %d, order in weak set: %s" % (magic(K, (z^43)^4), n43 in W))
cnt = Counter(magic(K, K.random_element()) for _ in range(2000))
print("random B, 2000 samples, magic number histogram:", dict(sorted(cnt.items())))
