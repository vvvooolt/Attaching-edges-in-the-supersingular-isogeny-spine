"""
v7b_ground.sage -- Stage 7b GROUND TRUTH (no class groups, no ideals).

For a prime l and each prime p of the input list, and each odd-t family D' of that l (logs/v7_odd_families.txt):
  * the roots of H_D' in F_p that are supersingular (Sage is_supersingular);
  * new l-edges: for each such root r, the roots r' != r of Phi_l(r, Y) mod p whose multiplicity exceeds the number of
    F_p-rational l-isogenies of E_r to curves with j = r' (Sage isogenies_prime_degree); the bridges of the family are
    the pairs {r, r'} with r' also a root of H_D'.  New neighbours outside the family are recorded as anomalies;
  * attachment of each bridge in the spine without new l-edges, with the j = 1728 gluing: components of F_p-classes
    under rational l-isogenies, twisting, and the gluing (both 1728 classes when p = 3 mod 4); k = j-graph distance;
  * levels of the endpoints (p = 3 mod 4): surface iff full rational 2-torsion.
Output per (p, family): n_bridges_obs, sorted multiset of 'A' / 'N<k>' tokens, levels, anomalies.
Usage: sage v7b_ground.sage L PRIMES.txt OUT.csv [FAMILYFILE]   (default logs/v7_odd_families.txt; the family
file's columns after l D' t f are ignored here)
"""
import sys, csv
from collections import deque
L = int(sys.argv[1])
FAMFILE = sys.argv[4] if len(sys.argv) > 4 else 'logs/v7_odd_families.txt'
FAMS = [tuple(int(x) for x in line.split()) for line in open(FAMFILE)
        if line.strip() and not line.startswith('#') and int(line.split()[0]) == L]
ps = [int(x) for x in open(sys.argv[2]).read().split()]
PHI = PolynomialRing(ZZ, 'x,y')(str(pari.polmodular(L)))
HS = {D: hilbert_class_polynomial(D) for (_, D, *_r) in FAMS}


def key(E):
    E = E.short_weierstrass_model()
    A, B = E.a4(), E.a6(); j = E.j_invariant()
    if A == 0:
        return (j, B.is_square())
    if B == 0:
        return (j, A.is_square())
    return (j, (3 * B / (2 * (1728 - j) * A)).is_square())


def level(E, p):
    return 'F' if p % 4 == 1 else ('S' if len(E.two_division_polynomial().roots()) == 3 else 'F')


rows = []
for p in ps:
    F = GF(p); R.<Y> = F[]
    PhiF = PHI.change_ring(F)
    ratcache = {}

    def rat(E):
        k = key(E)
        if k not in ratcache:
            ratcache[k] = [c.codomain() for c in E.isogenies_prime_degree(L)]
        return ratcache[k]

    comp_of, jadj = {}, {}

    def component(E):
        """j-vertex set of the spine component of j(E) (without new l-edges, with the 1728 gluing)."""
        k0 = key(E)
        if k0 in comp_of:
            return comp_of[k0]
        seen = {k0: E}; dq = deque([E]); glued = False; js = set()
        while dq:
            X = dq.popleft(); js.add(X.j_invariant())
            nb = rat(X)
            for Y2 in [X.quadratic_twist()] + nb:
                if key(Y2) not in seen:
                    seen[key(Y2)] = Y2; dq.append(Y2)
            for Y2 in nb:
                a, b = X.j_invariant(), Y2.j_invariant()
                if a != b:
                    jadj.setdefault(a, set()).add(b); jadj.setdefault(b, set()).add(a)
            if not glued and p % 4 == 3 and F(1728) in js:
                glued = True
                for Z in (EllipticCurve(F, [1, 0]), EllipticCurve(F, [-1, 0])):
                    if key(Z) not in seen:
                        seen[key(Z)] = Z; dq.append(Z)
        cid = frozenset(js)
        for k in seen:
            comp_of[k] = cid
        return cid

    def jdist(s, t):
        d = {s: 0}; dq = deque([s])
        while dq:
            u = dq.popleft()
            if u == t:
                return d[u]
            for v in jadj.get(u, ()):
                if v not in d:
                    d[v] = d[u] + 1; dq.append(v)
        return None

    for (l, D, t, f, *_rest) in FAMS:
        roots = [r for r, m in HS[D].change_ring(F).roots()]
        ss = [r for r in roots if r not in (F(0), F(1728)) and EllipticCurve(F, j=r).is_supersingular()]
        special = [r for r in roots if r in (F(0), F(1728))]
        bridges, anomalies, levels = set(), [], []
        for r in ss:
            E = EllipticCurve(F, j=r)
            mult = {rr: m for rr, m in R(PhiF(r, Y)).roots()}
            rc = {}
            for C in rat(E):
                rc[C.j_invariant()] = rc.get(C.j_invariant(), 0) + 1
            for rr, m in mult.items():
                extra = m - rc.get(rr, 0)
                if rr != r and extra > 0:
                    if rr in ss:
                        bridges.add(tuple(sorted((int(r), int(rr)))))
                    else:
                        anomalies.append('%d->%d x%d' % (r, rr, extra))
        toks = []
        for a, b in sorted(bridges):
            Ea, Eb = EllipticCurve(F, j=a), EllipticCurve(F, j=b)
            levels.append(level(Ea, p) + level(Eb, p))
            ca = component(Ea)
            if F(b) in ca:
                toks.append('N%d' % jdist(F(a), F(b)))
            else:
                toks.append('A')
        rows.append([p, p % 4, l, D, len(ss), len(bridges), ' '.join(sorted(toks)), ' '.join(sorted(set(levels))),
                     ' '.join(anomalies[:6]), ' '.join(str(s) for s in special)])
    sys.stdout.flush()
with open(sys.argv[3], 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['p', 'p_mod4', 'l', 'D', 'n_ss_roots', 'n_bridges_obs', 'bridges_obs', 'levels', 'anomalies', 'special_roots'])
    w.writerows(rows)
print("ground truth:", len(rows), "rows")
