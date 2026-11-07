"""lib.sage -- helpers shared by check_l2.sage and check_odd.sage.

Ground truth uses only curves over F_p and Sage's isogeny code (no class groups).
Predictions use only binary quadratic forms (PARI), never curves.
"""
from collections import deque

# ---------------------------------------------------------------- curves over F_p

def key(E):
    """Invariant of the F_p-isomorphism class of a supersingular E (j and twist class)."""
    E = E.short_weierstrass_model()
    A, B = E.a4(), E.a6()
    j = E.j_invariant()
    if A == 0:
        return (j, B.is_square())
    if B == 0:
        return (j, A.is_square())
    return (j, (3 * B / (2 * (1728 - j) * A)).is_square())


def models(F, j):
    """The two F_p-classes with invariant j (a curve and its quadratic twist)."""
    E = EllipticCurve(F, j=F(j))
    return [E, E.quadratic_twist()]


def rational_targets(E, l):
    """j-invariants of the targets of the F_p-rational l-isogenies of E (with multiplicity)."""
    return [phi.codomain().j_invariant() for phi in E.isogenies_prime_degree(l)]


def spine_minus_dist(F, j_start, j_goal, l):
    """Distance from j_start to j_goal in S^-_l: vertices are j-invariants in F_p, edges come from
    F_p-rational l-isogenies of either model.  Returns (distance or None, component size)."""
    dist = {j_start: 0}
    queue = deque([j_start])
    while queue:
        j = queue.popleft()
        for E in models(F, j):
            for j2 in rational_targets(E, l):
                if j2 not in dist:
                    dist[j2] = dist[j] + 1
                    queue.append(j2)
    return dist.get(j_goal), len(dist)


def level(E):
    """For p = 3 mod 4: 'surface' if all of E[2] is rational, else 'floor'."""
    return 'surface' if len(E(0).division_points(2)) == 4 else 'floor'


def phi_mult(F, l, j1, j2):
    """Multiplicity of j2 as a root of Phi_l(j1, Y) over F_p."""
    R = PolynomialRing(F, 'Y')
    Y = R.gen()
    P = R(classical_modular_polynomial(l)(F(j1), Y))
    return dict(P.roots()).get(F(j2), 0)


def rational_edge_exists(F, j1, j2, l):
    return any(F(t) == F(j2) for E in models(F, j1) for t in rational_targets(E, l))


# ---------------------------------------------------------------- class groups as forms

def one(D):
    return pari.qfbred(pari.qfbprimeform(D, 1))


def prime_form(D, q):
    """Class of a prime of norm q in the order of discriminant D (q must split)."""
    return pari.qfbred(pari.qfbprimeform(D, q))


def mul(a, b):
    return pari.qfbred(pari.qfbcomp(a, b))


def pw(a, n):
    return pari.qfbred(pari.qfbpow(a, n))


def cyclic(g, D):
    """Powers of g: dict str(form) -> exponent, and the order of g."""
    e, x, seen = 0, one(D), {}
    while str(x) not in seen:
        seen[str(x)] = e
        x, e = mul(x, g), e + 1
    return seen, e


def dlog(c, table):
    return table.get(str(c))


def nrm(e, o):
    """Distance of e to 0 in Z/o."""
    e %= o
    return min(e, o - e)


def half(e, o):
    """e/2 modulo the odd number o."""
    return (e * inverse_mod(2, o)) % o


# ---------------------------------------------------------------- class-group action on the floor

def _split_y(P, R, f):
    A, B = R(0), R(0)
    if P in P.parent().base_ring():
        return R(P), R(0)
    for (ex, ey), c in P.dict().items():
        t = c * R.gen()**ex * f**(ey // 2)
        if ey % 2:
            B += t
        else:
            A += t
    return A, B


def frobenius_eigenvalue(E, K, q, p):
    """Eigenvalue of Frobenius (mod q) on the rational kernel with kernel polynomial K.
    Computed in F_p[x]/(K): x^p = x([lam]P) and y^p/y = y([lam]P)/y."""
    R = K.parent()
    Q = R.quotient(K)
    f = R.gen()**3 + E.a4() * R.gen() + E.a6()
    xp = Q(R.gen())**p
    yratio = Q(f)**((p - 1) // 2)
    for lam in range(1, q):
        mx, my = E.multiplication_by_m(lam)
        nA, _ = _split_y(mx.numerator(), R, f)
        dA, _ = _split_y(mx.denominator(), R, f)
        if Q(dA) == 0 or xp != Q(nA) / Q(dA):
            continue
        _, nB = _split_y(my.numerator(), R, f)
        dA, _ = _split_y(my.denominator(), R, f)
        if yratio == Q(nB) / Q(dA):
            return lam
    return None


def act(E, q, lam, p):
    """[l]E for l = (q, pi - lam): the rational q-isogeny on whose kernel Frobenius acts by lam."""
    E = E.short_weierstrass_model()
    hits = [phi.codomain() for phi in E.isogenies_prime_degree(q)
            if frobenius_eigenvalue(E, phi.kernel_polynomial(), q, p) == lam % q]
    assert len(hits) == 1
    return hits[0]


def label_floor(p):
    """Label every floor curve (p = 3 mod 4) by its class in Cl(Z[sqrt(-p)]), starting from
    E0: y^2 = x^3 + x.  Returns dict key(E) -> (form, j), the class number, and the number of
    inconsistencies (a curve reached twice with different labels)."""
    F, D = GF(p), -4 * p
    h = int(pari(D).qfbclassno())
    e = one(D)
    gens, sub = [], {str(e)}
    for q in prime_range(3, 400):
        if len(sub) == h:
            break
        if kronecker(-p, q) != 1:
            continue
        lam = int(GF(q)(-p).sqrt())
        for l in (lam, q - lam):
            gens.append((q, l, pari.qfbred(pari.Qfb(q, 2 * l, (l * l + p) // q))))
        sub, todo = {str(e)}, [e]
        while todo:
            g = todo.pop()
            for _, _, f in gens:
                n = mul(g, f)
                if str(n) not in sub:
                    sub.add(str(n))
                    todo.append(n)
    E0 = EllipticCurve(F, [1, 0])
    label = {key(E0): (e, E0.j_invariant())}
    frontier, bad = [E0], 0
    while frontier:
        nxt = []
        for E in frontier:
            lab = label[key(E)][0]
            for q, l, f in gens:
                E2 = act(E, q, l, p)
                k2, new = key(E2), mul(lab, f)
                if k2 in label:
                    bad += label[k2][0] != new
                else:
                    label[k2] = (new, E2.j_invariant())
                    nxt.append(E2)
        frontier = nxt
    return label, h, bad
