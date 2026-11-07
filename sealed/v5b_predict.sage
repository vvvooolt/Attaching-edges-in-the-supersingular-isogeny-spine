"""
v5b_predict.sage -- Phase 5b PREDICTIONS (class-group data only), to be hashed before any
verification run.  For l = 3 and the family D' = -35 (H_-35 roots in F_p and supersingular:
(5/p) = 1, (-7/p) = -1, p not in {2,3,5,7,23}).

Uses only binary quadratic forms of discriminant -4p (the order Z[sqrt(-p)]; the endpoints
are floor curves by Prop. 12.x(c) / Lemma 12.3):
  H3 = <prime form above 3> if p = 2 mod 3, else trivial;  o3 = |H3|;
  attach_pred  = [q5] not in H3 and [q7] not in H3                       (Theorem A3)
  k_pred       = min ||e_q|| over q in {5,7} with [q_q] = [q_3]^{e_q} in H3 (Corollary C3)
  forced_pred  = 7 * 3^floor(o3/2) < p                                    (Corollary C3)
  for p = 3 mod 4 and p < POS_LIMIT: the unordered pair of floor classes of the edge,
  {c+^((h-1)/2), c-^((h-1)/2)} with c+- = [q5 q7^{+-1}], each up to inversion (Theorem B3(a)),
  written as reduced forms (a,b,c) with b >= 0 chosen canonically.
Usage: sage v5b_predict.sage PRIMES.txt OUT.csv
"""
import sys, csv
POS_LIMIT = 30000


def canon(f):
    """Canonical representative of {f, f^-1} (reduced forms)."""
    g = pari.qfbred(pari.Qfb(f[0], -f[1], f[2]))
    a, b = (int(f[0]), int(f[1]), int(f[2])), (int(g[0]), int(g[1]), int(g[2]))
    return min(a, b)


def predict(p):
    D = -4 * p
    red = lambda f: pari.qfbred(f)
    mul = lambda u, v: pari.qfbred(pari.qfbcomp(u, v))
    inv = lambda f: pari.qfbred(pari.Qfb(f[0], -f[1], f[2]))
    one = pari.Qfb(1, 0, p)
    h = int(pari.qfbclassno(D))
    q5, q7 = red(pari.qfbprimeform(D, 5)), red(pari.qfbprimeform(D, 7))
    Hd = {str(one): 0}
    if p % 3 == 2:
        g = red(pari.qfbprimeform(D, 3)); x, k = g, 1
        while str(x) not in Hd:
            Hd[str(x)] = k; x = mul(x, g); k += 1
    o3 = len(Hd)
    rim = lambda e: min(e % o3, (-e) % o3)
    in5, in7 = str(q5) in Hd, str(q7) in Hd
    att = int(not in5 and not in7)
    ks = [rim(Hd[str(f)]) for f, b in ((q5, in5), (q7, in7)) if b]
    k = min(ks) if ks else ''
    forced = int(7 * 3 ** (o3 // 2) < p)
    pos = ''
    if p % 4 == 3 and p < POS_LIMIT:
        assert h % 2 == 1
        cp = pari.qfbpow(mul(q5, q7), (h - 1) // 2)
        cm = pari.qfbpow(mul(q5, inv(q7)), (h - 1) // 2)
        pos = ' '.join(sorted('(%d,%d,%d)' % canon(red(f)) for f in (cp, cm)))
    return [p, p % 3, p % 4, h, o3, att, int(in5), int(in7), k, forced, pos]


if len(sys.argv) >= 3 and sys.argv[1].endswith('.txt'):
    ps = [int(x) for x in open(sys.argv[1]).read().split()]
    with open(sys.argv[2], 'w') as fh:
        w = csv.writer(fh)
        w.writerow(['p', 'p_mod3', 'p_mod4', 'h_4p', 'o3', 'attach_pred', 'q5_in_H3', 'q7_in_H3',
                    'k_pred', 'forced_attach_pred', 'edge_floor_classes_pred'])
        for p in ps:
            w.writerow(predict(p))
