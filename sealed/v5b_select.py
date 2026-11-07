"""v5b_select.py -- prime lists for Phase 5b (written before any verification run).
  logs/v5b_inrange_primes.txt : ALL family primes p < 10^6 ((5/p)=1, (-7/p)=-1, p not in {2,3,5,7,23})
  logs/v5b_oos_primes.txt     : 200 family primes in (10^6, 4*10^6), seed 5202609, 20 equal-width
                                bands, 10 per band."""
import random
import common as C
fam = lambda p: p not in (2, 3, 5, 7, 23) and C.leg(5, p) == 1 and C.leg(-7, p) == -1
inr = [p for p in range(5, 10**6) if C.isprime(p) and fam(p)]
open('logs/v5b_inrange_primes.txt', 'w').write(' '.join(map(str, inr)) + '\n')
random.seed(5202609)
LO, HI, NB = 10**6, 4 * 10**6, 20
w = (HI - LO) // NB
oos = []
for b in range(NB):
    a = LO + b * w
    ps = [p for p in range(a, a + w) if fam(p) and C.isprime(p)]
    oos += sorted(random.sample(ps, 10))
open('logs/v5b_oos_primes.txt', 'w').write(' '.join(map(str, oos)) + '\n')
print(len(inr), "in-range family primes;", sum(1 for p in inr if p % 3 == 2), "with p = 2 mod 3;",
      len(oos), "OOS primes", min(oos), max(oos))
