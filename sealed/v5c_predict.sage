"""
v5c_predict.sage -- Phase 5c PREDICTIONS (class-group data only), to be sealed before any
verification run.  l = 3, family D' = -20: H_-20 = x^2 - 1264000 x - 681472000 has two supersingular
roots in F_p iff p = 11, 19 (mod 20)  (so p = 3 mod 4 and both volcano levels exist).

Theorem A_20 (proof.pdf Section 13): one endpoint is on the floor, the other on the surface; the edge is
NOT attaching iff p = 2 (mod 3) and [l5] lies in H3f = <[l3]> in Cl(Z[sqrt(-p)]) (disc -4p).
Predicted per prime, using only binary quadratic forms of discriminants -4p (floor) and -p (surface):
  o3f, o3s   = orders of the prime forms above 3 (1 if p = 1 mod 3);
  in_f, in_s = [l5] in <[l3]> on the floor / on the surface (in_f implies in_s);
  attach_pred = not in_f;
  k_pred     = ||e5f/2||_{o3f} + ||e5s/2||_{o3s}  (distance through j = 1728), when not attaching,
               where [l5] = [l3]^e5 and /2 is inversion of 2 modulo the odd order;
  forced_attach_pred = 5 * 3^floor(o3f/2) < p;
  levels_pred = 'SF' (one endpoint on each level).
Usage: sage v5c_predict.sage PRIMES.txt OUT.csv
"""
import sys, csv


def sub_and_exp(D, q, target):
    """Cyclic subgroup <prime form above q> (as dict str(form) -> exponent) and the exponent of target."""
    red = lambda f: pari.qfbred(f)
    one = red(pari.qfbpow(pari.qfbprimeform(D, q), 0))
    g = red(pari.qfbprimeform(D, q))
    Hd, x, k = {str(one): 0}, g, 1
    while str(x) not in Hd:
        Hd[str(x)] = k; x = red(pari.qfbcomp(x, g)); k += 1
    return len(Hd), Hd.get(str(red(target)))


def rim(e, o):
    return min(e % o, (-e) % o)


def predict(p):
    Df, Ds = -4 * p, -p
    hf, hs = int(pari.qfbclassno(Df)), int(pari.qfbclassno(Ds))
    if p % 3 == 1:
        return [p, p % 3, p % 8, hf, hs, 1, 1, 0, 0, 1, '', 1, 'SF']
    o3f, e5f = sub_and_exp(Df, 3, pari.qfbprimeform(Df, 5))
    o3s, e5s = sub_and_exp(Ds, 3, pari.qfbprimeform(Ds, 5))
    in_f, in_s = e5f is not None, e5s is not None
    att = int(not in_f)
    k = ''
    if not att:
        half = lambda e, o: e * ((o + 1) // 2)
        k = rim(half(e5f, o3f), o3f) + rim(half(e5s, o3s), o3s)
    forced = int(5 * 3 ** (o3f // 2) < p)
    return [p, p % 3, p % 8, hf, hs, o3f, o3s, int(in_f), int(in_s), att, k, forced, 'SF']


if len(sys.argv) >= 3 and sys.argv[1].endswith('.txt'):
    ps = [int(x) for x in open(sys.argv[1]).read().split()]
    with open(sys.argv[2], 'w') as fh:
        w = csv.writer(fh)
        w.writerow(['p', 'p_mod3', 'p_mod8', 'h_4p', 'h_p', 'o3f', 'o3s', 'q5_in_H3f', 'q5_in_H3s',
                    'attach_pred', 'k_pred', 'forced_attach_pred', 'levels_pred'])
        for p in ps:
            w.writerow(predict(p))
