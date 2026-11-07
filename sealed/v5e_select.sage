"""v5e_select.sage -- out-of-sample primes for Theorem A_gen (Phase 5e), chosen before any computation on them:
for each prime-degree family (l, D') with l in {5, 7}: 50 family primes in (10^5, 10^6), seed 5202611,
10 per band, 5 equal-width bands.  Family: kronecker(D', p) = -1 and H_D' has two distinct roots mod p.
Writes logs/v5e_oos_primes_<l>_<|D'|>.txt."""
import random
random.seed(int(5202611))
for l, D in [(5, -51), (5, -91), (7, -115), (7, -187)]:
    H = hilbert_class_polynomial(D)
    out = []
    for b in range(5):
        a, c = 10**5 + b * 180000, 10**5 + (b + 1) * 180000
        band = [int(p) for p in prime_range(a, c) if kronecker(D, p) == -1
                and len(H.change_ring(GF(p)).roots()) == 2]
        out += sorted(random.sample(band, 10))
    open('logs/v5e_oos_primes_%d_%d.txt' % (l, -D), 'w').write(' '.join(map(str, out)) + '\n')
    print(l, D, len(out), out[:3], out[-2:])
