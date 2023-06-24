from sage.all import *  #  to avoid circular imports triggered by pyhdme in old versions of sage
from pyhdme import (
    modular_igusa_from_igusa_clebsch,
    igusa_clebsch_from_modular_igusa,
    siegel_modeq_isog_invariants_Q_wrapper,
    siegel_modeq_2step_isog_invariants_Q_wrapper,
)
from sage.all import (
    GCD,
    GF,
    Graph,
    HyperellipticCurve,
    HyperellipticCurve_from_invariants,
    Integers,
    PolynomialRing,
    Primes,
    QQ,
    ZZ,
    cached_function,
    kronecker_symbol,
    matrix,
    magma,
    prime_range,
    prod,
    sqrt,
    squarefree_part,
    vector,
)
import sys

from subprocess import Popen, PIPE

def rescale(c, I, weights):
    return vector([c**i * j for i, j in zip(weights, I)])
def make_integral(I, weights):
    return rescale(LCM([elt.denominator() for elt in I]), I, weights)

def modular_invariants(C):
    return modular_igusa_from_igusa_clebsch(
        make_integral(C.change_ring(QQ).igusa_clebsch_invariants(), (1,2,3,5))
    )




@cached_function
def Lpolynomial(C, p):
    try:
        Cp = C.change_ring(GF(p))
    except ValueError:
        return None
    if Cp.genus() == C.genus():
        return Cp.frobenius_polynomial()
    else:
        return None


@cached_function
def discriminant(C):
    f, h = C.hyperelliptic_polynomials()
    return (h**2 + 4 * f).discriminant()

@cached_function
def quadratic_twist(C, d):
    if d == 1:
        return C
    f, h = C.hyperelliptic_polynomials()
    return HyperellipticCurve(d*(h**2 + 4 * f))



def HyperellipticCurve_from_modular_invariants(minv, reduced=True):
    ic = igusa_clebsch_from_modular_igusa(minv)
    try:
        C = HyperellipticCurve_from_invariants(ic, reduced=False)
    except ZeroDivisionError:
        C = magma(ic).ChangeUniverse(QQ).HyperellipticCurveFromIgusaClebsch().sage()
    C = ReducedMinimalWeierstrassModel(C)
    # lazy way to normalize the igusa
    newminv = modular_invariants(C)

    # find the scalar d that relates minv and newminv in P(4,6,10,12)
    weights = [4, 6, 10, 12]
    coordinates = [i for i, elt in enumerate(minv) if elt != 0]
    assert coordinates == [i for i, elt in enumerate(newminv) if elt != 0]
    assert len(coordinates) > 0
    c = coordinates[0]
    d = (newminv[c]/minv[c]).nth_root(weights[c])
    if reduced:
        # we don't care about spoiling the model over Q, and thus we can call ReducedWamelenModel
        return ReducedMinimalWeierstrassModel(ReducedWamelenModel(quadratic_twist(C, d)))
    else:
        return quadratic_twist(C, d)

def possible_isogenous_quadratic_twists(C, bad_primes, Lpolynomial_origin, bound=2000):
    C = C.change_ring(Integers())
    ZZx = PolynomialRing(ZZ, "x")
    x = ZZx.gen()
    bad_primes += [-1]
    bad_primes += ZZ(discriminant(C)).prime_divisors()
    bad_primes = sorted(set(bad_primes))
    twistdata = []
    splittingdata = []
    if bound is None:
        primes = Primes()
    else:
        primes = prime_range(bound)

    for p in primes:
        if p in bad_primes:
            continue
        efc = ZZx(Lpolynomial(C, p))
        efo = ZZx(Lpolynomial_origin(p))
        efctwist = efc.subs({x: -x})
        if efo not in [efc, efctwist]:
            return []
        if efc == efctwist:
            continue
        if efc == efo:
            twistdata.append(0)
        else:
            twistdata.append(1)

        splittingdata.append(
            [1 if kronecker_symbol(q, p) == -1 else 0 for q in bad_primes]
        )
        if matrix(GF(2), splittingdata).rank() == len(bad_primes):
            break
    M = matrix(GF(2), splittingdata).transpose().stack(vector(GF(2), twistdata))
    k = M.nrows()
    solutions = [v[: k - 1] for v in M.kernel().basis() if v[k - 1] == 1]
    r = [
        prod(bad_primes[i] for i, e in enumerate(sol) if e != 0) for sol in solutions
    ]
    return r


def isogenous_curves(C, invariants, reduced=True, known_models={}):
    invariants = list(map(tuple, invariants))
    # returns curves in the same order
    bad_primes = ZZ(discriminant(C)).prime_divisors()
    Lpolynomial_origin = lambda p: Lpolynomial(C, p)
    Cinv = [
        known_models.get(tuple(elt), HyperellipticCurve_from_modular_invariants(elt, reduced=reduced).change_ring(ZZ))
        for elt in invariants
    ]
    twists = [
        possible_isogenous_quadratic_twists(elt, bad_primes, Lpolynomial_origin)
        for elt in Cinv
    ]
    assert all(len(elt) == 1 for elt in twists)
    Cnonred = [quadratic_twist(c, t[0]) for c, t in zip(Cinv, twists)]
    if reduced:
        res = []
        for inv, nonredmodel in zip(invariants, Cnonred):
            res.append(known_models.get(inv, ReducedMinimalWeierstrassModel(nonredmodel)))
    else:
        res = Cnonred
    assert [modular_invariants(elt) for elt in res] == invariants
    return res


def gpwrapper(C, command):
    f, h = C.hyperelliptic_polynomials()
    sfh = [f.list(), h.list()]
    pipe = Popen(["/usr/local/bin/gp","-q"], stdin=PIPE, stdout=PIPE, encoding="utf8")
    cmd = "{print(iferr(%s(apply(Polrev,%s)), E, \"None\"))}" % (command, sfh);
    out, _ = pipe.communicate(input=f"{cmd};quit()")
    if out == "None\n":
        return C
    R = PolynomialRing(Integers(), 'x')
    S = f.parent()
    newf, newh = map(S, map(R, out.strip('[]\n').split(",")))
    return HyperellipticCurve(newf, newh)

def MinimalWeierstrassModel(C):
    return gpwrapper(C, 'hyperellminimalmodel')

def ReducedModel(C):
    return gpwrapper(C, 'hyperellred')


def ReducedMinimalWeierstrassModel(C):
    # minimize discriminant
    C0 = MinimalWeierstrassModel(C)
    # minimize coefficients
    C1 = ReducedModel(C0)
    return C1

def ReducedWamelenModel(C):
    """
    Given a hyperelliptic curve C over Q, returns a reduced and partially minimized model of some quadratic twist of C
    """
    return magma(C).ReducedWamelenModel().sage()

def IsIsomorphic(C0, C1):
    return magma(C0).IsIsomorphic(magma(C1)).sage()


def isogeny_graph_invariants(ics, ells, verbose=0, threads=1):
    G = Graph(weighted=True)
    queue = {tuple(ics[0])}
    done = set()
    steps = [
        siegel_modeq_isog_invariants_Q_wrapper,
        siegel_modeq_2step_isog_invariants_Q_wrapper,
    ]
    while queue:
        elt = queue.pop()
        G.add_vertex(elt)
        if elt in done:
            continue
        if verbose:
            print(f"handling = {elt}")
        for ell in ells:
            if verbose:
                print(f"ell = {ell}")
            for i, s in enumerate(steps):
                for elt1 in s(elt, ell, verbose=verbose, threads=threads):
                    elt1 = tuple(elt1)
                    G.add_vertex(elt1)
                    G.add_edge(elt, elt1, ell ** (i + 1))
                    if verbose:
                        print(f"adding {elt1}")
                    queue.add(elt1)
        if verbose:
            print(f"done with = {elt}")
        done.add(elt)
        if verbose:
            print(f"queue = {queue}")
            print(f"done = {done}")
    assert set(map(tuple, ics)).issubset(G.vertices(sort=False))
    return G


def isogeny_graph(C, conductor, ells=None, verbose=0, threads=1, reduced=True, known_models={}):
    if ells is None:
        ells = reducible_ell(C, conductor)
        if verbose:
            print(f"ells = {ells}")
    m_inv = modular_invariants(C)
    G = isogeny_graph_invariants([m_inv], ells, verbose=verbose, threads=threads)
    invs = G.vertices(sort=True)
    curves = isogenous_curves(C, invs, reduced=reduced, known_models=known_models)
    invs_to_curves = dict([(i, c) for i, c in zip(invs, curves)])
    G.relabel(invs_to_curves)
    return G, ells


def isogeny_graph_from_line(line, verbose=0, threads=1, reduced=True):
    cond, c, Lhash, input_ics, input_eqns = line.split(":")
    cond = int(cond)
    input_ics, input_eqns = map(eval, [input_ics, input_eqns])
    ZZx = PolynomialRing(ZZ, 'x')
    known_models = dict(zip(map(tuple, input_ics), [HyperellipticCurve(ZZx(f), ZZx(h)) for f, h in input_eqns]))
    R = PolynomialRing(ZZ, "x")
    # pick the first curve
    C = HyperellipticCurve(*map(R, input_eqns[0]))
    G, ells = isogeny_graph(C, cond, verbose=verbose, threads=threads, reduced=reduced, known_models=known_models)
    curves = G.vertices(key=modular_invariants)
    eqns = [[pol.list() for pol in elt.hyperelliptic_polynomials()] for elt in curves]
    invs = list(map(modular_invariants, curves))
    assert invs == sorted(invs)
    assert set(map(tuple, input_ics)).issubset(set(invs))
    invs = list(map(list, invs))
    # to be sure we have everything in the same order
    M = G.weighted_adjacency_matrix(vertices=curves)
    return (
        cond,
        c,
        Lhash,
        input_ics,
        input_eqns,
        ells,
        invs,
        eqns,
        list(map(list, M.rows())),
    )


def reducible_ell(C, conductor, bound=500):
    bad_primes = ZZ(discriminant(C)).prime_divisors()
    sqfree_cond = squarefree_part(conductor)
    d = sqrt(conductor / sqfree_cond)
    Zd = Integers(d)

    def onedim_subrep(Lp, p, f):
        x = Lp.parent().gen()
        return Lp.resultant(x**f - 1)

    def twodim_subrep1(Lp, p, f):
        x = Lp.parent().gen()
        ap = -Lp[3]
        bp = Lp[2]
        Q1p = (bp * x - 1 - p**2 * x**2) * (p * x + 1) ** 2 - ap**2 * p * x**2
        return Q1p.resultant(x**f - 1)

    def twodim_subrep2(Lp, p, f):
        x = Lp.parent().gen()
        ap = -Lp[3]
        bp = Lp[2]
        Q2p = (bp * x - p - p * x**2) * (x + 1) ** 2 - ap**2 * x**2
        return Q2p.resultant(x**f - 1)

    tests = [onedim_subrep, twodim_subrep1, twodim_subrep2]
    res = [0 for _ in tests]
    for p in prime_range(bound):
        if p in bad_primes:
            continue
        f = Zd(p).multiplicative_order()
        Lp = Lpolynomial(C, p)
        if not Lp:
            continue
        for i, t in enumerate(tests):
            if res[i] != 1:
                res[i] = GCD(res[i], t(Lp, p, f))
        if all(r == 1 for r in res):
            break
    ells = set(sum([list(elt.prime_divisors()) for elt in res], []))
    if not ells:
        return []
    polys = [Lpolynomial(C, p) for p in prime_range(bound) if p not in bad_primes]
    # there are some good primes that we also loose due to bad choice of model
    polys = [elt for elt in polys if elt]
    new_ells = []
    for ell in ells:
        polys_ell = [elt.change_ring(GF(ell)) for elt in polys]
        if not any(
            elt.change_ring(GF(ell)).is_irreducible()
            for elt in polys_ell
            if elt[0] % ell != 0
        ):
            new_ells.append(ell)
    return new_ells

def is_interactive():
    import __main__ as main
    return not hasattr(main, '__file__')

if __name__ == "__main__" and not is_interactive():
    if not len(sys.argv) >= 2:
        print(r"""
    USAGE:
         sage -python {cmd} <input line> <threads=1> <verbose=0> <reduced=0>

         where the last 3 positional arguments are optional
         the <input line> takes the following format
            "conductor:<optional>:<optional>:[list of modular invariants]:[list of pairs of coefficients]"
        for example
            "704:704.a:2159575534599616393:[[-1856,129536,-44,-2948]]:[[[-1,-1,-2,1,1,-2],[1,1,1,1]]]"

    EXAMPLE
        """.format(cmd=sys.argv[0]) )
        sys.exit()
    threads = 1 if len(sys.argv) <= 2 else int(sys.argv[2])
    verbose = 0 if len(sys.argv) <= 3 else int(sys.argv[3])
    reduced = True if len(sys.argv) <= 4 else bool(int(sys.argv[4]))
    out = isogeny_graph_from_line(sys.argv[1], verbose=verbose, threads=threads, reduced=reduced)
    print(":".join(map(str, out)).replace(" ", ""))
