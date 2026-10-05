# Candidate field sizes for a randomized (per-namespace) curve:
# odd m in [85, 127], rejecting composite m (GHS descent to a subfield), and for
# prime m the lowest irreducible trinomial/pentanomial plus the cost of finding
# a random curve  y^2 + xy = x^3 + x^2 + B  with a prime-order subgroup of
# index two (#E = 2r).
#
#   nix develop .#math -c sage sage/field_sizes.sage
#
import time
R.<x> = GF(2)[]
def lowest_poly(m):
    for k in range(1, m // 2 + 1):
        if (x^m + x^k + 1).is_irreducible(): return x^m + x^k + 1
    for k3 in range(3, m):
        for k2 in range(2, k3):
            for k1 in range(1, k2):
                f = x^m + x^k3 + x^k2 + x^k1 + 1
                if f.is_irreducible(): return f
print(" m  bytes  factor   reduction poly               log2 r  counts  ms/count")
for m in map(ZZ, range(85, 129, 2)):
    if not m.is_prime():
        p = factor(m)[0][0]
        print("%3d  %2d    %-8s -- descent to GF(2^%d), n = %d%s" % (m, ceil(m / 8), factor(m), m // p, p,
              " (magic <= 3 for every b)" if p == 3 else ""))
        continue
    f = lowest_poly(m)
    K.<z> = GF(2^m, modulus=f)
    t0 = time.time(); c = 0
    while True:
        B = K.random_element(); c += 1
        n = EllipticCurve(K, [1, 1, 0, 0, B]).order()
        if n % 4 == 2 and (n // 2).is_pseudoprime(): break
    print("%3d  %2d    prime    %-28s %.1f   %3d     %.1f" % (m, ceil(m / 8), str(f), float(log(n // 2, 2)), c, 1000 * (time.time() - t0) / c))
