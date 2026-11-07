"""check_odd.sage -- the odd-l families of Section 6 (and through them Theorems D and E).

Usage:  sage check_odd.sage LO HI [FAMILY ...]      FAMILY in: 3,-35 5,-51 5,-91 7,-115 7,-187 3,-20 3,-32
                                                    (default: all seven)

For each family (l, D') and each prime p in [LO, HI) above the family's bound at which H_D' has two
distinct supersingular roots j1, j2 in F_p:

Ground truth (curves only):
  * j1 j2 is a double edge of the l-isogeny graph and no model of j1 has a rational l-isogeny to j2;
  * the levels of the endpoints (p = 3 mod 4);
  * attachment and length k = d(j1, j2) in S^-_l, by breadth-first search (vertices are
    j-invariants, so the two j = 1728 classes are glued automatically).
Prediction (forms only):
  * prime degrees (Cor 6.1): both ends on the floor; attaching iff [l_{n+}], [l_{n-}] not in
    <[l_l]> in Cl(Z[sqrt -p]); length min ||x|| over [l_q] = [l_l]^x; bound q l^k >= p (Thm E);
  * D' = -20 (Cor 6.4): one end on each level; not attaching iff p = 2 mod 3 and [l_5] in <[l_3]>;
    k = ||e/2||_{o_f} + ||e/2||_{o_s}; 5 * 3^||e|| >= p;
  * D' = -32 (Cor 6.6): one end on each level; not attaching iff p = 2 mod 3 and [p_2] in <[l_3]>
    in Cl(O_K); k = ||3e/2||_o + ||e/2||_o.
"""
import sys
load('lib.sage')

LO, HI = int(sys.argv[1]), int(sys.argv[2])
#            (l, D')     : (n+, n-) for prime degrees, bound on p (Remark 5.4, Table 2)
FAMILIES = {(3, -35): ((7, 5), 315), (5, -51): ((17, 3), 1275), (5, -91): ((13, 7), 2275),
            (7, -115): ((23, 5), 5635), (7, -187): ((17, 11), 9163),
            (3, -20): (None, 180), (3, -32): (None, 288)}
chosen = [tuple(int(v) for v in a.split(',')) for a in sys.argv[3:]] or list(FAMILIES)


def subgroup(D, l, p):
    """<[l_l]> in the class group of discriminant D (trivial if l does not split)."""
    g = prime_form(D, l) if kronecker(-p, l) == 1 else one(D)
    return cyclic(g, D)


grand = []
for (l, Dp) in chosen:
    degs, bound = FAMILIES[(l, Dp)]
    HD = hilbert_class_polynomial(Dp)
    n = n_att = 0
    fails = []
    for p in prime_range(max(LO, bound + 1), HI):
        F = GF(p)
        roots = sorted(r for r, _ in HD.change_ring(F).roots())
        if len(roots) != 2 or not EllipticCurve(F, j=roots[0]).is_supersingular():
            continue
        j1, j2 = roots
        n += 1
        bad = []
        if phi_mult(F, l, j1, j2) != 2 or rational_edge_exists(F, j1, j2, l):
            bad.append('edge not a new double edge')
        d, _ = spine_minus_dist(F, j1, j2, l)
        attaches = d is None
        n_att += attaches
        lev = sorted(level(EllipticCurve(F, j=j)) for j in (j1, j2)) if p % 4 == 3 else None
        Df = -4 * p

        if degs is not None:                                   # Cor 6.1
            if lev is not None and lev != ['floor', 'floor']:
                bad.append('levels %s' % lev)
            T, o = subgroup(Df, l, p)
            xs = [dlog(prime_form(Df, q), T) for q in degs]
            in_H = [(q, nrm(x, o)) for q, x in zip(degs, xs) if x is not None]
            if (not in_H) != attaches:
                bad.append('attachment')
            elif not attaches:
                q, k = min(in_H, key=lambda t: t[1])
                if k != d:
                    bad.append('length %d != %d' % (k, d))
                if q * l**d < p:
                    bad.append('Theorem E bound')

        elif Dp == -20:                                        # Cor 6.4
            if lev != ['floor', 'surface']:
                bad.append('levels %s' % lev)
            Tf, of = subgroup(Df, 3, p)
            e = dlog(prime_form(Df, 5), Tf)
            pred_att = not (p % 3 == 2 and e is not None)
            if pred_att != attaches:
                bad.append('attachment')
            elif not attaches:
                Ts, os_ = subgroup(-p, 3, p)
                k = nrm(half(e, of), of) + nrm(half(e, os_), os_)
                if k != d:
                    bad.append('length %d != %d' % (k, d))
                if 5 * 3**nrm(e, of) < p:
                    bad.append('bound 5*3^m >= p')

        elif Dp == -32:                                        # Cor 6.6
            if lev != ['floor', 'surface']:
                bad.append('levels %s' % lev)
            Ts, o = subgroup(-p, 3, p)
            e = dlog(prime_form(-p, 2), Ts)
            pred_att = not (p % 3 == 2 and e is not None)
            if pred_att != attaches:
                bad.append('attachment')
            elif not attaches:
                k = nrm(half(3 * e, o), o) + nrm(half(e, o), o)
                if k != d:
                    bad.append('length %d != %d' % (k, d))

        if bad:
            fails.append((p, bad))
        print('l=%d D=%d p=%d attaches=%s k=%s %s' % (l, Dp, p, attaches, d,
              'OK' if not bad else 'FAIL ' + ', '.join(bad)), flush=True)
    grand.append((l, Dp, n, n_att, fails))

print()
for l, Dp, n, n_att, fails in grand:
    print("SUMMARY l=%d D'=%d: %d primes in [%d, %d), %d attaching; failures: %s"
          % (l, Dp, n, LO, HI, n_att, fails if fails else 'none'))
