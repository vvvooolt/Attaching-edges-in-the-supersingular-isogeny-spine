"""v5d_predict.sage -- Phase 5d PREDICTIONS from Conjectures A32 / C32 (logs/v5d_conjecture.txt), class-group
data only (binary quadratic forms of discriminant -p).  Sealed before any verification on these primes.
Usage: sage v5d_predict.sage PRIMES.txt OUT.csv"""
import sys, csv


def predict(p):
    D = -p
    h = int(pari.qfbclassno(D))
    if p % 3 == 1:
        return [p, p % 3, h, 1, '', 1, 1, 'SF']
    g = pari.qfbprimeform(D, 3); H = {}; x = pari.qfbred(pari.qfbpow(g, 0)); k = 0
    while str(x) not in H:
        H[str(x)] = k; x = pari.qfbred(pari.qfbcomp(x, g)); k += 1
    o = len(H); e = H.get(str(pari.qfbred(pari.qfbprimeform(D, 2))))
    n = lambda y: min(y % o, (-y) % o)
    half = (o + 1) // 2
    att = int(e is None)
    kp = '' if att else n(3 * e * half) + n(e * half)
    forced = int(8 * 3 ** (o // 2) < p + 1)
    return [p, p % 3, h, o, kp if not att else '', att, forced, 'SF']


if len(sys.argv) >= 3 and sys.argv[1].endswith('.txt'):
    ps = [int(x) for x in open(sys.argv[1]).read().split()]
    with open(sys.argv[2], 'w') as fh:
        w = csv.writer(fh)
        w.writerow(['p', 'p_mod3', 'h_p', 'o3', 'k_pred', 'attach_pred', 'forced_attach_pred', 'levels_pred'])
        for p in ps:
            w.writerow(predict(p))
