"""
v7b_predict.sage -- Stage 7b PREDICTIONS from Theorem A_odd (proof.pdf Part II), class-group data only
(binary quadratic forms of discriminant -4p; no curves, no modular polynomials).

Families (l, D', t, f, e+, e-) are read from logs/v7_odd_families.txt.  For each prime p and each family of that l:
  * if some prime q | e+ e- has kronecker(-4p, q) != 1: 0 bridges predicted;
  * else: bridges <-> pairs ({a, a-bar}, {b, b-bar}) of primitive ideals of norms e+, e- (all sign choices, up to
    conjugation of each factor); for each: attach = [a] not in H_l and [b] not in H_l, H_l = <prime form above l>
    (trivial if l does not split); k = min rim distance of those of [a], [b] in H_l.
  Output per (p, family): n_bridges and the sorted multiset of (attach, k) as 'A' / 'N<k>' tokens.
Usage: sage v7b_predict.sage L PRIMES.txt OUT.csv
"""
import sys, csv, itertools
L = int(sys.argv[1])
FAMS = [tuple(int(x) for x in line.split()) for line in open('logs/v7_odd_families.txt')
        if line.strip() and not line.startswith('#') and int(line.split()[0]) == L]
ps = [int(x) for x in open(sys.argv[2]).read().split()]


def red(f):
    return pari.qfbred(f)


def ideal_classes(e, Dp):
    """Classes of primitive ideals of norm e (odd), one per conjugation pair; None if some prime does not split."""
    fac = list(factor(e))
    if any(kronecker(Dp, q) != 1 for q, _ in fac):
        return None
    base = [red(pari.qfbpow(pari.qfbprimeform(Dp, q), k)) for q, k in fac]
    out = []
    for signs in itertools.product([1, -1], repeat=len(fac) - 1):
        g = base[0]
        for s, b in zip(signs, base[1:]):
            g = red(pari.qfbcomp(g, b if s == 1 else pari.Qfb(b[0], -b[1], b[2])))
        out.append(g)
    return out


rows = []
for p in ps:
    Dp = -4 * p
    one = pari.Qfb(1, 0, p)
    H = {str(one): 0}
    if kronecker(Dp, L) == 1:
        g = red(pari.qfbprimeform(Dp, L)); x, k = g, 1
        while str(x) not in H:
            H[str(x)] = k; x = red(pari.qfbcomp(x, g)); k += 1
    o = len(H)
    rim = lambda e: min(e % o, (-e) % o)
    for (l, D, t, f, ep, em) in FAMS:
        A, B = ideal_classes(ep, Dp), ideal_classes(em, Dp)
        toks = []
        if A is not None and B is not None:
            for a in A:
                for b in B:
                    ea, eb = H.get(str(red(a))), H.get(str(red(b)))
                    if ea is None and eb is None:
                        toks.append('A')
                    else:
                        toks.append('N%d' % min(rim(e) for e in (ea, eb) if e is not None))
        rows.append([p, p % 4, l, D, t, f, ep, em, o, len(toks), ' '.join(sorted(toks))])
with open(sys.argv[3], 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['p', 'p_mod4', 'l', 'D', 't', 'f', 'e_plus', 'e_minus', 'o_l', 'n_bridges_pred', 'bridges_pred'])
    w.writerows(rows)
print("predictions:", len(rows), "rows")
