# Survey of GLS-style curves over GF(2^(2m)) = GF(2^m)[u]/(u^2+u+1), m odd,
# near one 64-bit word per base-field element. Rejects m with a small odd
# factor p (GHS descent to GF(2^(2m/p)), n = p); for prime m reports the
# first sparse b giving #E' = 2*prime.
#
#   nix develop .#math -c sage sage/gls_survey.sage
#
# GLS-style curves E'/GF(q^2), q = 2^m, m odd:  y^2 + xy = x^3 + u x^2 + b,
# u^2 + u + 1 = 0, b in GF(q).  #E'(GF(q^2)) = (q-1)^2 + t^2, t = trace of
# y^2 + xy = x^3 + b over GF(q).
R.<x> = GF(2)[]

def lowest_poly(m):
    for k in range(1, m // 2 + 1):
        if (x^m + x^k + 1).is_irreducible():
            return "x^%d+x^%d+1" % (m, k), x^m + x^k + 1
    for k3 in range(3, m):
        for k2 in range(2, k3):
            for k1 in range(1, k2):
                f = x^m + x^k3 + x^k2 + x^k1 + 1
                if f.is_irreducible():
                    return "x^%d+x^%d+x^%d+x^%d+1" % (m, k3, k2, k1), f

def first_prime_order(m, f):
    q = 2^m
    K.<w> = GF(q, modulus=f)
    cands = [(w^k, "w^%d" % k) for k in range(1, m)] + [(1 + w^k, "1+w^%d" % k) for k in range(1, m)] + [(1 + w^j + w^k, "1+w^%d+w^%d" % (j, k)) for k in range(2, m) for j in range(1, k)]
    for tries, (b, name) in enumerate(cands, 1):
        t = q + 1 - EllipticCurve(K, [1, 0, 0, 0, b]).order()
        n = (q - 1)^2 + t^2
        if n % 2 == 0 and (n // 2).is_pseudoprime():
            return K, b, name, n, tries
    return K, None, None, None, len(cands)

print("m   factor(m)   reduction poly               first b    log2 r   rho(neg)  rho(neg+psi)")
for m in map(ZZ, range(41, 65, 2)):
    name, f = lowest_poly(m)
    fm = factor(m)
    if not m.is_prime():
        print("%-3d %-11s %-28s  -- subfield GF(2^%d) of GF(2^%d): descent with n = %d"
              % (m, fm, name, 2 * m // fm[0][0], 2 * m, fm[0][0]))
        continue
    K, b, bname, n, tries = first_prime_order(m, f)
    r = n // 2
    print("%-3d %-11s %-28s  %-9s  %.2f   2^%.1f    2^%.1f"
          % (m, fm, name, bname, float(log(r, 2)),
             float(log(sqrt(pi.n() * r / 4), 2)), float(log(sqrt(pi.n() * r / 8), 2))))
    if m == 61:
        K2 = GF(2^(2 * m), 'v')
        emb = K.hom([f.change_ring(K2).roots()[0][0]])
        u = (x^2 + x + 1).change_ring(K2).roots()[0][0]
        assert EllipticCurve(K2, [1, u, 0, 0, emb(b)]).order() == n
        print("    m=61: direct point count over GF(2^122) agrees with (q-1)^2 + t^2")
