"""v5e_predict.sage -- Phase 5e PREDICTIONS from Theorem A_gen (proof.pdf Thm 15.2), class-group data only
(binary quadratic forms of discriminant -4p).  Sealed before any ground truth on the out-of-sample primes.
For family (l, D', n+, n-):  H_l = <prime form above l> if (-p/l) = 1 else trivial;
attach_pred = neither prime form above n+- lies in H_l;  k_pred = min rim distance of those in H_l;  o_pred = |H_l|.
Usage: sage v5e_predict.sage L D N1,N2 PRIMES.txt OUT.csv"""
import sys, csv
l, D = int(sys.argv[1]), int(sys.argv[2]); DEGS = [int(x) for x in sys.argv[3].split(',')]
ps = [int(x) for x in open(sys.argv[4]).read().split()]
with open(sys.argv[5], 'w') as fh:
    w = csv.writer(fh); w.writerow(['p', 'l', 'D', 'o_pred', 'n_in_H', 'attach_pred', 'k_pred'])
    for p in ps:
        Dp = -4 * p; one = pari.Qfb(1, 0, p); H = {str(one): 0}
        if kronecker(-p, l) == 1:
            g = pari.qfbred(pari.qfbprimeform(Dp, l)); x, k = g, 1
            while str(x) not in H:
                H[str(x)] = k; x = pari.qfbred(pari.qfbcomp(x, g)); k += 1
        o = len(H); ins = ''; ks = []
        for n in DEGS:
            e = H.get(str(pari.qfbred(pari.qfbprimeform(Dp, n))))
            ins += '1' if e is not None else '0'
            if e is not None:
                ks.append(min(e % o, (-e) % o))
        w.writerow([p, l, D, o, ins, int(ins == '00'), min(ks) if ks else ''])
