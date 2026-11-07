"""v7_oos_select.sage -- Phase 7 out-of-sample primes (chosen before any computation on them), seed 7202609.
For l = 3: 40 primes in (10^6, 3*10^6) (all l = 3 families were already checked below 10^6);
for l = 5, 7: 40 primes in (10^5, 5*10^5);  l = 11: (58564, 3*10^5);  l = 13: (114244, 3*10^5)  (all above 4 l^4 and above
every in-sample range).  Uniform random sample of the primes in the range.  Writes logs/v7_oos/primes_<l>.txt."""
import random
random.seed(int(7202609))
RANGES = {3: (10**6, 3 * 10**6), 5: (10**5, 5 * 10**5), 7: (10**5, 5 * 10**5), 11: (58565, 3 * 10**5), 13: (114245, 3 * 10**5)}
for l in (3, 5, 7, 11, 13):
    lo, hi = RANGES[l]
    ps = [int(p) for p in prime_range(lo, hi)]
    S = sorted(random.sample(ps, 40))
    open('logs/v7_oos/primes_%d.txt' % l, 'w').write(' '.join(map(str, S)) + '\n')
    print(l, S[:3], S[-2:], sum(1 for p in S if p % 4 == 3), 'with p = 3 mod 4')
