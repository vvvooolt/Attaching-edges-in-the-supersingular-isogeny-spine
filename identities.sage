"""identities.sage -- every fact marked "computed" in the paper, plus the numbers in its examples.

Usage:  sage identities.sage            (prints PASS/FAIL per item and a summary)
"""
load('lib.sage')
R.<X> = ZZ[]
results = []


def check(name, ok):
    results.append((name, bool(ok)))
    print('%-4s %s' % ('PASS' if ok else 'FAIL', name), flush=True)


H = lambda D: hilbert_class_polynomial(D)(X)
Phi = lambda l: classical_modular_polynomial(l)

# ---- Lemma 3.1: H_-15
H15 = H(-15)
check('H_-15 = X^2 + 191025 X - 121287375', H15 == X**2 + 191025*X - 121287375)
check('disc H_-15 = 3^6 5^3 7^4 13^2', H15.discriminant() == 3**6 * 5**3 * 7**4 * 13**2)
check('H_-15(0) = -3^6 5^3 11^3', H15(0) == -3**6 * 5**3 * 11**3)
check('H_-15(1728) = 3^6 7^4 11^2', H15(1728) == 3**6 * 7**4 * 11**2)

# ---- Table 1 and Remark 3.12: Res_Y(Phi_l, dPhi_l/dY) as products of class polynomials
S.<x, y> = ZZ[]
D5 = [-99, -96, -91, -84, -64, -51, -36, -24]
D7 = [-195, -192, -187, -180, -171, -160, -132, -115, -96, -75, -52, -48, -40, -20]
TABLE1 = {2: (-4, {-15: 2, -7: 2, -4: 1, -3: 2}),
          3: (-27, {-35: 2, -32: 2, -20: 2, -11: 2, -8: 2, -4: 2, -3: 2}),
          5: (-3125, dict([(D, 2) for D in D5] + [(-19, 2), (-16, 2), (-11, 2), (-4, 4), (-3, 4)])),
          7: (-823543, dict([(D, 2) for D in D7] + [(-27, 2), (-24, 2), (-19, 2), (-12, 2), (-4, 4), (-3, 6)]))}
for l, (c, mult) in TABLE1.items():
    P = Phi(l)(x, y)
    res = R(P.resultant(P.derivative(y), y).univariate_polynomial().list())
    expected = c * prod(H(D)**m for D, m in mult.items())
    check('Table 1, l = %d: Res_l = %d * product of H_D' % (l, c), res == expected)
check('Remark 3.12: Res_2 = -4 X^2 (X-1728) (X+3375)^2 H_-15^2',
      TABLE1[2][0] * prod(H(D)**m for D, m in TABLE1[2][1].items())
      == -4 * X**2 * (X - 1728) * (X + 3375)**2 * H15**2)

# ---- Table 2 and Lemma 4.5 hypotheses: families with prime degrees
for l, Dp, t0, npm, bound in [(3, -35, 1, (7, 5), 315), (5, -51, 7, (17, 3), 1275), (5, -91, 3, (13, 7), 2275),
                              (7, -115, 9, (23, 5), 5635), (7, -187, 3, (17, 11), 9163)]:
    ok = (Dp == t0**2 - 4 * l**2 and Dp.is_fundamental_discriminant() and QuadraticField(Dp).class_number() == 2
          and (2 * l + t0, 2 * l - t0) == npm and all(is_prime(n) for n in npm)
          and (1 - Dp) // 4 > l and abs(Dp) > l**2 and t0 % l != 0
          and bound == max(l**2 * abs(Dp), 4 * l**2))
    check('Table 2, l = %d, D\' = %d' % (l, Dp), ok)

# ---- Section 6.2: the three families for l = 3
for Dp, discH, H0, H1728, t0 in [(-35, 2**32 * 5**3 * 7**2 * 23**2, -2**30 * 5**3, 2**12 * 7**2 * 19**2 * 31**2, 1),
                                 (-20, 2**18 * 5**3 * 13**2 * 17**2, -2**12 * 5**3 * 11**3, -2**16 * 11**2 * 19**2, 4),
                                 (-32, 2**9 * 5**6 * 7**4 * 13**2 * 29**2, 2**6 * 5**6 * 23**3, -2**6 * 7**4 * 23**2 * 31**2, 2)]:
    HD = H(Dp)
    principal = BinaryQF([1, Dp % 2, (Dp % 2 - Dp) // 4])      # the principal form of discriminant D'
    principal_rep_3 = principal.solve_integer(3) is not None
    check('l = 3, D\' = %d: disc, H(0), H(1728), h = 2, t0, prime above 3 not principal' % Dp,
          HD.discriminant() == discH and HD(0) == H0 and HD(1728) == H1728 and HD.degree() == 2
          and Dp == t0**2 - 36 and not principal_rep_3)
check('Res(H_-32, H_-8) = -2^6 5^6 7^4 13^2', H(-32).resultant(H(-8)) == -2**6 * 5**6 * 7**4 * 13**2)
check('root fields Q(sqrt 5), Q(sqrt 5), Q(sqrt 2)',
      [squarefree_part(H(D).discriminant()) for D in (-35, -20, -32)] == [5, 5, 2])

# ---- Corollary 6.3: the H_-140 identity
P2 = Phi(2)(x, y)
r = R(H(-35)(x).resultant(P2, x).univariate_polynomial().list())
check('Res_X(H_-35(X), Phi_2(X,Y)) = H_-140(Y)', r == H(-140))
check('Res(H_-35, H_-140) = 2^48 5^21 7^14 19^8 23^4 31^4 37^4 41^2 43^4 53^2 59^2 61^3 67^2',
      H(-35).resultant(H(-140)) == 2**48 * 5**21 * 7**14 * 19**8 * 23**4 * 31**4 * 37**4 * 41**2 * 43**4
      * 53**2 * 59**2 * 61**3 * 67**2)

# ---- Example 1.1 and Section 3.8
first = [p for p in prime_range(60, 400) if p % 120 in (71, 119)]
check('first five primes 71,119 mod 120 have <[p_2]> = Cl(O_K)',
      first[:5] == [71, 191, 239, 311, 359]
      and all(cyclic(prime_form(-p, 2), -p)[1] == int(pari(-p).qfbclassno()) for p in first[:5]))
T, o = cyclic(prime_form(-431, 2), -431)
check('p = 431: h = 21, ord [p_2] = 7, [q_3], [q_5] not in <[p_2]>, roots {150, 189}',
      int(pari(-431).qfbclassno()) == 21 and o == 7 and dlog(prime_form(-431, 3), T) is None
      and dlog(prime_form(-431, 5), T) is None
      and sorted(r for r, _ in H15.change_ring(GF(431)).roots()) == [150, 189])
T, o = cyclic(prime_form(-6551, 2), -6551)
check('p = 6551: h = 117, ord [p_2] = 39, [q_3] = [p_2]^(+-15), [q_5] not in <[p_2]>, (479^2 + 25 p)/4 = 3 * 2^15',
      int(pari(-6551).qfbclassno()) == 117 and o == 39 and nrm(dlog(prime_form(-6551, 3), T), o) == 15
      and dlog(prime_form(-6551, 5), T) is None and (479**2 + 25 * 6551) // 4 == 3 * 2**15)
p = 327599
T, o = cyclic(prime_form(-p, 2), -p)
h = int(pari(-p).qfbclassno())
q3, q5 = prime_form(-p, 3), prime_form(-p, 5)
xs = [nrm(dlog(pw(c, (h - 1) // 2), T), o) for c in (mul(q3, q5), mul(q3, pari.qfbpow(q5, -1)))]
check('p = 327599: h = ord [p_2] = 963, ||e_3|| = 354, ||e_5|| = 14, positions {170, 184}, (81 + p)/4 = 5 * 2^14',
      h == o == 963 and nrm(dlog(q3, T), o) == 354 and nrm(dlog(q5, T), o) == 14 and sorted(xs) == [170, 184]
      and (81 + p) // 4 == 5 * 2**14 and p == 20 * 2**14 - 81)
check('sharpness: 7 * 3^6 = 5103 (p = 5099), 17 * 7^4 = 40817 (p = 40801)', 7 * 3**6 == 5103 and 17 * 7**4 == 40817)

print('\nSUMMARY identities: %d checks, failures: %s'
      % (len(results), [n for n, ok in results if not ok] or 'none'))
