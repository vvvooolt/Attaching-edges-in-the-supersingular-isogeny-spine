"""v5c_select.py -- prime list for Phase 5c (written before any verification run):
  logs/v5c_primes.txt : ALL family primes p < 10^6 for D' = -20, i.e. p = 11, 19 (mod 20)."""
import common as C
ps = [p for p in range(11, 10**6) if p % 20 in (11, 19) and C.isprime(p)]
open('logs/v5c_primes.txt', 'w').write(' '.join(map(str, ps)) + '\n')
print(len(ps), "family primes;", sum(1 for p in ps if p % 3 == 2), "with p = 2 mod 3")
