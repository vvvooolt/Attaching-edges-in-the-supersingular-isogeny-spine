"""
v7_predict_all.sage -- Stage 7f PREDICTIONS for every family, p = 3 (mod 4), from class-group data only.

Abstract model (no curves, no modular polynomials):
  * floor curves  F(b), b in Cl_f = Cl(Z[sqrt -p])  (forms of disc -4p), F(b) = [b] E0f, E0f: y^2 = x^3 + x;
    surface curves S(c), c in Cl_s = Cl(O_K)  (forms of disc -p),   S(c) = [c] E0s, E0s: y^2 = x^3 - x;
    twist: b -> b^-1 (both levels: the bases are their own twists);
    im: Cl_f -> Cl_s, [a] -> [a O_K] (a group homomorphism, built on generators: the odd prime form
    (q, 2 lam, .) of disc -4p maps to the prime form (q, B, .) of disc -p with B = lam mod q, B odd).
  * rational steps: odd prime q: multiply by the prime form above q (either orientation), same level;
    2-moves: ascend F(b) -> S(im b); descend S(c) -> F(b) for each b with im b = c;
    horizontal (p = 7 mod 8 only) S(c) -> S(c * P2^{+-1}), P2 = (2, 1, (p+1)/8).
    A walk is cyclic iff it never backtracks (P2 then P2^-1; ascend then descend to the same floor class;
    descend then ascend).
  * R_e(X) = endpoints of cyclic rational walks of degree e from X (odd part: primitive ideals of norm e_odd,
    all orientations; 2-part: all cyclic 2-walks of length v2(e)).
Family (l, D', t, f, n+, e+, n-, e-):  m0 = |D'| if D' = 1 mod 4 else |D'|/4.
  roots  = { X : iota X in R_{m0}(X) }, restricted by the optimal-order rule:
           D' = 1 mod 4 -> X on the floor;   D' = 0 mod 4 -> not (X on the floor and m0 = 3 mod 4);
  bridge partner of X = the root M (j(M) != j(X)) with  M in R_{e+}(X), iota M in R_{e-}(X)  or with + and - swapped;
  attachment (proved rules): same level: not attaching iff pos(M) in H pos(X)^{+-1}; different levels: not attaching iff
  pos(X) in H and pos(M) in H (H = <prime form above l> on that level, trivial if l does not split);
  k = min rim distance (same level) or sum of the two rim distances to 1728 (different levels).
Output per (p, family): n_roots_pred (as j-invariants), n_bridges_pred, sorted multiset 'A' / 'N<k>', flags.
Usage: sage v7_predict_all.sage L PRIMES.txt OUT.csv FAMILYFILE
"""
import sys, csv, itertools
L = int(sys.argv[1])
ps = [int(x) for x in open(sys.argv[2]).read().split()]
FAMS = []
for line in open(sys.argv[4]):
    if line.strip() and not line.startswith('#'):
        v = [int(x) for x in line.split()]
        if v[0] != L:
            continue
        if len(v) == 6:          # odd-family file: l D t f e+ e-  (n+- computed as needed)
            l, D, t, f, ep, em = v
        else:
            l, D, t, f, npl, ep, nmi, em = v
        FAMS.append((l, D, t, f, ep, em))


def red(f):
    return pari.qfbred(f)


def key(f):
    return (int(f[0]), int(f[1]), int(f[2]))


def predict(p):
    Df, Ds = -4 * p, -p
    onef = pari.Qfb(1, 0, p)
    ones = pari.Qfb(1, 1, (1 + p) // 4)
    hf, hs = int(pari.qfbclassno(Df)), int(pari.qfbclassno(Ds))
    mulf = lambda a, b: red(pari.qfbcomp(a, b))
    inv = lambda a: red(pari.Qfb(a[0], -a[1], a[2]))
    # generators of Cl_f by small split odd primes, with images in Cl_s
    gens = []
    sub = {key(onef): ones}
    for q in prime_range(3, 2000):
        if len(sub) == hf:
            break
        if kronecker(-p, q) != 1:
            continue
        lam = int(GF(q)(-p).sqrt())
        gf = red(pari.Qfb(q, 2 * lam, (lam * lam + p) // q))
        B = lam if lam % 2 == 1 else lam + q
        gs = red(pari.Qfb(q, B, (B * B + p) // (4 * q)))
        gens.append((gf, gs))
        sub = {key(onef): ones}; todo = [(onef, ones)]
        while todo:
            a, b = todo.pop()
            for gf_, gs_ in gens:
                n = mulf(a, gf_)
                if key(n) not in sub:
                    sub[key(n)] = mulf(b, gs_); todo.append((n, sub[key(n)]))
    assert len(sub) == hf
    im = sub                                    # key(floor form) -> surface form
    fiber = {}
    for kf, s in im.items():
        fiber.setdefault(key(s), []).append(kf)
    assert len(fiber) == hs
    floor_forms = {kf: None for kf in im}
    P2 = red(pari.Qfb(2, 1, (1 + p) // 8)) if p % 8 == 7 else None

    def form(kk):
        return pari.Qfb(*kk)

    # odd ideal classes of norm e (all orientations), on each level
    def odd_ideals(e, level):
        fac = list(factor(e)) if e > 1 else []
        D = Df if level == 'F' else Ds
        if any(kronecker(D, q) != 1 for q, _ in fac):
            return None
        opts = []
        for q, k in fac:
            g = red(pari.qfbprimeform(D, q))
            gk = red(pari.qfbpow(g, k))
            opts.append([gk, inv(gk)])
        one = onef if level == 'F' else ones
        out = []
        for ch in itertools.product(*opts):
            x = one
            for y in ch:
                x = red(pari.qfbcomp(x, y))
            out.append(key(x))
        return out

    def walks2(X, a):
        """endpoints of cyclic 2-walks of length a from X = (level, key)."""
        res = []

        def rec(cur, prev_move, k):
            if k == 0:
                res.append(cur); return
            lev, c = cur
            if lev == 'F':
                if prev_move == 'down':
                    return                      # ascending right after descending backtracks
                rec(('S', key(im[c])), 'up', k - 1)
            else:
                # descend
                for kf in fiber[c]:
                    if prev_move == 'up' and kf == prev_from:
                        continue
                    rec2(('F', kf), 'down', k - 1)
                if P2 is not None:
                    for sgn in (1, -1):
                        if prev_move == ('h', -sgn):
                            continue
                        g = P2 if sgn == 1 else inv(P2)
                        rec(('S', key(red(pari.qfbcomp(form(c), g)))), ('h', sgn), k - 1)
        # prev_from tracking for the 'up then down' backtrack
        def rec2(cur, mv, k):
            rec(cur, mv, k)
        prev_from = None
        if X[0] == 'F' and a > 0:
            # first step from floor is the ascent; remember the start to forbid descending back
            lev, c = X
            prev_from = c
            rec(('S', key(im[c])), 'up', a - 1)
            # rec uses prev_from for the descent check
            return res
        rec(X, None, a)
        return res

    def R(X, e):
        eo, a = e, 0
        while eo % 2 == 0:
            eo //= 2; a += 1
        ids = odd_ideals(eo, X[0])
        if ids is None:
            return set()
        out = set()
        for idk in ids:
            Y = (X[0], key(red(pari.qfbcomp(form(X[1]), form(idk)))))
            for Z in walks2(Y, a):
                out.add(Z)
        return out

    def twist(X):
        return (X[0], key(inv(form(X[1]))))

    # H and exponents on each level
    def cyc(level):
        D = Df if level == 'F' else Ds
        one = onef if level == 'F' else ones
        Hd = {key(one): 0}
        if kronecker(D, L) == 1:
            g = red(pari.qfbprimeform(D, L)); x, k = g, 1
            while key(x) not in Hd:
                Hd[key(x)] = k; x = red(pari.qfbcomp(x, g)); k += 1
        return Hd
    HF, HS = cyc('F'), cyc('S')
    allX = [('F', kf) for kf in im] + [('S', ks) for ks in fiber]
    rows = []
    for (l, D, t, f, ep, em) in FAMS:
        Dodd = (D % 4 == 1)
        m0 = -D if Dodd else -D // 4
        roots = []
        for X in allX:
            if Dodd and X[0] != 'F':
                continue
            if (not Dodd) and X[0] == 'F' and m0 % 4 == 3:
                continue
            if twist(X) in R(X, m0):
                roots.append(X)
        # j-classes: X ~ twist(X)
        jcls = {}
        for X in roots:
            jcls.setdefault(frozenset([X, twist(X)]), X)
        flags = []
        base_hits = [X for X in roots if X[1] in (key(onef), key(ones))]
        if base_hits:
            flags.append('root_at_1728')
        partner = {}
        for J, X in jcls.items():
            A, B = R(X, ep), R(X, em)
            cands = set()
            for J2, M in jcls.items():
                if J2 == J:
                    continue
                if (M in A and twist(M) in B) or (twist(M) in A and M in B) or \
                   (M in B and twist(M) in A) or (twist(M) in B and M in A):
                    cands.add(J2)
            partner[J] = cands
        bridges = set()
        for J, c in partner.items():
            if len(c) != 1:
                flags.append('partners=%d' % len(c))
            for J2 in c:
                bridges.add(frozenset([J, J2]))
        toks = []
        for Bg in bridges:
            J1, J2 = list(Bg)
            X, M = jcls[J1], jcls[J2]
            HX = HF if X[0] == 'F' else HS
            HM = HF if M[0] == 'F' else HS
            if X[0] == M[0]:
                Hd = HX; o = len(Hd)
                q1 = key(red(pari.qfbcomp(form(M[1]), inv(form(X[1])))))
                q2 = key(red(pari.qfbcomp(form(M[1]), form(X[1]))))
                es = [Hd[q] for q in (q1, q2) if q in Hd]
                toks.append('A' if not es else 'N%d' % min(min(e % o, (-e) % o) for e in es))
            else:
                if X[1] in HX and M[1] in HM:
                    rx = min(HX[X[1]] % len(HX), (-HX[X[1]]) % len(HX))
                    rm = min(HM[M[1]] % len(HM), (-HM[M[1]]) % len(HM))
                    toks.append('N%d' % (rx + rm))
                else:
                    toks.append('A')
        levels = ' '.join(sorted(set(''.join(sorted(jcls[J][0] for J in Bg)) for Bg in bridges)))
        rows.append([p, p % 8, l, D, t, f, ep, em, m0, len(jcls), len(bridges), ' '.join(sorted(toks)), levels,
                     ' '.join(sorted(set(flags)))])
    return rows


def predict_1mod4(p):
    # p = 1 (mod 4): one level, Cl = Cl(Z[sqrt -p]) of even order.  CONJECTURAL enumeration (Conj. A_1mod4):
    # bridges <-> pairs ({a+, conj}, {a-, conj}) of primitive odd ideals of norms e+-_odd, with 2-parts t^{c+-},
    # c+- = v2(e+-) (0 bridges if some c > 1), existing iff all odd primes split and the twisting class
    # [a+ a- t^{c+ + c-}] is a non-square (twisting classes are never squares for p = 1 mod 4);
    # attach iff [a+ t^{c+}] and [a- t^{c-}] are not in H = <prime form above l>.
    D = -4 * p
    one = pari.Qfb(1, 0, p)
    h = int(pari.qfbclassno(D))
    Tt = red(pari.Qfb(2, 2, (1 + p) // 2))
    # all classes and squares
    gens = []; cls = {key(one): one}
    for q in prime_range(2, 3000):
        if len(cls) == h:
            break
        if q == 2:
            g = Tt
        elif kronecker(D, q) == 1:
            g = red(pari.qfbprimeform(D, q))
        else:
            continue
        gens.append(g)
        todo = list(cls.values())
        while todo:
            a = todo.pop()
            for g_ in gens:
                n = red(pari.qfbcomp(a, g_))
                if key(n) not in cls:
                    cls[key(n)] = n; todo.append(n)
    assert len(cls) == h
    squares = {key(red(pari.qfbcomp(a, a))) for a in cls.values()}
    Hd = {key(one): 0}
    if kronecker(D, L) == 1:
        g = red(pari.qfbprimeform(D, L)); x, k = g, 1
        while key(x) not in Hd:
            Hd[key(x)] = k; x = red(pari.qfbcomp(x, g)); k += 1
    o = len(Hd)
    rows = []
    for (l, D_, t, f, ep, em) in FAMS:
        def split(e):
            c = 0
            while e % 2 == 0:
                e //= 2; c += 1
            return e, c
        (oa, ca), (ob, cb) = split(ep), split(em)
        toks = []
        ok = ca <= 1 and cb <= 1
        def pats(e):
            fac = list(factor(e)) if e > 1 else []
            if any(kronecker(D, q) != 1 for q, _ in fac):
                return None
            base = [red(pari.qfbpow(pari.qfbprimeform(D, q), k)) for q, k in fac]
            if not base:
                return [one]
            out = []
            for signs in itertools.product([1, -1], repeat=len(base) - 1):
                x = base[0]
                for sg, b in zip(signs, base[1:]):
                    x = red(pari.qfbcomp(x, b if sg == 1 else red(pari.Qfb(b[0], -b[1], b[2]))))
                out.append(x)
            return out
        A, B = (pats(oa), pats(ob)) if ok else (None, None)
        if A is not None and B is not None:
            tw = lambda x, c: red(pari.qfbcomp(x, Tt)) if c else x
            c0 = key(red(pari.qfbcomp(pari.qfbcomp(A[0], B[0]), pari.qfbpow(Tt, ca + cb))))
            if c0 not in squares:
                m0 = -D_ if D_ % 4 == 1 else -D_ // 4
                mu = 1 if (D_ % 4 == 1 or m0 % 4 == 3) else 2
                rep = mu if (oa > 1 and ob > 1) else mu // 2
                for a in A:
                    for b in B:
                      for _rep in range(rep):
                        ga, gb = key(tw(a, ca)), key(tw(b, cb))
                        es = [Hd[x] for x in (ga, gb) if x in Hd]
                        toks.append('A' if not es else 'N%d' % min(min(e % o, (-e) % o) for e in es))
    # (rep: multiplicity mu = number of curves of each t-pair {X, [t]X} that are roots -- 1 for odd D' or m0 = 3 mod 4,
    #  2 otherwise; a pattern gives mu bridges if both e+ and e- have odd prime factors, mu/2 if one of them has none.)
        rows.append([p, p % 8, l, D_, t, f, ep, em, '', '', len(toks), ' '.join(sorted(toks)), 'FF' if toks else '',
                     'conj_1mod4'])
    return rows


out = []
for p in ps:
    if p < 11:
        continue
    out += predict(p) if p % 4 == 3 else predict_1mod4(p)
with open(sys.argv[3], 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['p', 'p_mod8', 'l', 'D', 't', 'f', 'e_plus', 'e_minus', 'm0', 'n_roots_pred', 'n_bridges_pred',
                'bridges_pred', 'levels_pred', 'flags'])
    w.writerows(out)
print("predictions:", len(out))
