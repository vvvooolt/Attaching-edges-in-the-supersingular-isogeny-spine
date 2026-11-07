"""check_l2.sage -- Theorems A, B, C (l = 2) against graphs computed from scratch.

Usage:  sage check_l2.sage LO HI [BMAX]

For every prime 60 < p < HI, p >= LO, with p = 11, 14 mod 15 (the primes where H_-15 has two
supersingular roots j1, j2 in F_p):

Ground truth (curves only, no class groups):
  * the edge j1 j2 is a double edge of the 2-isogeny graph and no F_p-model of j1 has a rational
    2-isogeny to j2 (it is new; Cor 3.6);
  * a model of j1 has F_p-rational isogenies of degrees 3 and 5 to models of j2 (Lemma 3.7);
  * for p = 3 mod 4 both endpoints are on the floor (Cor 3.6);
  * attachment and length k = d(j1, j2) - 2 by breadth-first search in S^-_2.
Prediction (binary quadratic forms only):
  * Theorem A: attaching iff [q3], [q5] are not in H_2;
  * Theorem C / 3.17: the length k, and the bound 5 * 2^k >= (p+1)/4; Cor 3.18;
  * Theorem B (p = 3 mod 4, p < BMAX, default 5000): every floor curve is labelled by its class in
    Cl(Z[sqrt -p]) by walking with Frobenius-eigenvalue isogenies from y^2 = x^3 + x, and
    {j1, j2} must equal {j(c+^((h'-1)/2) E0), j(c-^((h'-1)/2) E0)}.
Prints one line per prime and a summary; any failure is listed.
"""
import sys
load('lib.sage')

LO, HI = int(sys.argv[1]), int(sys.argv[2])
BMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 5000
H15 = [1, 191025, -121287375]
fails, n, n_att, n_B = [], 0, 0, 0

for p in prime_range(max(LO, 61), HI):
    if p % 15 not in (11, 14):
        continue
    n += 1
    F = GF(p)
    R = PolynomialRing(F, 'X')
    X = R.gen()
    j1, j2 = sorted(r for r, _ in (X**2 + 191025 * X - 121287375).roots())
    bad = []

    # ---- ground truth
    if not EllipticCurve(F, j=j1).is_supersingular():
        bad.append('not supersingular')
    if phi_mult(F, 2, j1, j2) != 2 or rational_edge_exists(F, j1, j2, 2):
        bad.append('edge not a new double edge')
    if not (rational_edge_exists(F, j1, j2, 3) and rational_edge_exists(F, j1, j2, 5)):
        bad.append('no rational 3- and 5-isogenies')
    if p % 4 == 3 and not (level(EllipticCurve(F, j=j1)) == level(EllipticCurve(F, j=j2)) == 'floor'):
        bad.append('endpoint not on the floor')
    d, size = spine_minus_dist(F, j1, j2, 2)
    attaches = d is None
    k_obs = None if attaches else d - 2

    # ---- prediction from forms
    DK = -p if p % 4 == 3 else -4 * p
    q3, q5 = prime_form(DK, 3), prime_form(DK, 5)
    if p % 8 == 7:
        g = prime_form(DK, 2)
    elif p % 4 == 1:
        g = prime_form(DK, 2)          # the class of t = (2, 1 + sqrt -p), of order 2
    else:
        g = one(DK)                    # 2 is inert: H_2 = 1
    T, o = cyclic(g, DK)
    e3, e5 = dlog(q3, T), dlog(q5, T)
    pred_att = e3 is None and e5 is None
    if pred_att != attaches:
        bad.append('Theorem A')
    else:
        n_att += attaches

    if not attaches:
        if p % 8 == 7:
            if e3 is not None and e5 is not None:          # folded component
                h = int(pari(-p).qfbclassno())
                x = dlog(pw(mul(q3, q5), (h - 1) // 2), T)
                y = dlog(pw(mul(q3, pari.qfbpow(q5, -1)), (h - 1) // 2), T)
                k_pred = abs(nrm(x, o) - nrm(y, o))
            else:                                          # stacked component
                k_pred = nrm(e3 if e3 is not None else e5, o)
            if k_pred != k_obs:
                bad.append('length %s != %s' % (k_pred, k_obs))
        if 5 * 2**k_obs < (p + 1) / 4:
            bad.append('Theorem C bound')
    if p % 8 == 7 and 5 * 2**((o - 1) // 2) < (p + 1) / 4 and not attaches:
        bad.append('Cor 3.18')

    # ---- Theorem B by labelling the floor
    if p % 4 == 3 and p < BMAX:
        n_B += 1
        label, hf, incons = label_floor(p)
        jof = {str(f): j for f, j in label.values()}
        Df = -4 * p
        c_plus = mul(prime_form(Df, 3), prime_form(Df, 5))
        c_minus = mul(prime_form(Df, 3), pari.qfbpow(prime_form(Df, 5), -1))
        pred = {jof.get(str(pw(c, (hf - 1) // 2))) for c in (c_plus, c_minus)}
        if incons or len(label) != hf or pred != {j1, j2}:
            bad.append('Theorem B')

    if bad:
        fails.append((p, bad))
    print('p=%d attaches=%s k=%s %s' % (p, attaches, k_obs, 'OK' if not bad else 'FAIL ' + ', '.join(bad)),
          flush=True)

print('SUMMARY l=2: %d primes in [%d, %d), %d attaching, Theorem B checked at %d primes; failures: %s'
      % (n, LO, HI, n_att, n_B, fails if fails else 'none'))
