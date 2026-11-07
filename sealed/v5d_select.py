"""v5d_select.py -- fresh prime list for Phase 5d: ALL family primes p = 7 (mod 8) with 10^5 < p < 10^6
(disjoint from the exploratory range p < 10^5)."""
import common as C
ps = [p for p in range(10**5 + 1, 10**6) if p % 8 == 7 and C.isprime(p)]
open('logs/v5d_fresh_primes.txt', 'w').write(' '.join(map(str, ps)) + '\n')
print(len(ps), "fresh family primes")
