"""Corpus generator. assemble_inputs() adds fixture/differential entries to
corpus_curves.json; compute_oracles() (Task 3) adds `expect` blocks and emits
tests/test_isogenyprimes.m and tests/run_differential.m.

Usage:
    sage -python tests/generate_corpus.py            # full run: writes the
        corpus with `expect` blocks and tests/test_isogenyprimes.m +
        tests/run_differential.m (minutes: the oE oracle screens each (curve,
        ell) to a certain-negative before any construction, and compute_oracles
        forks across entries; the orchestrator runs this).
    sage -python tests/generate_corpus.py --smoke    # deterministic subset,
        emits /tmp/smoke_*.m only, never touches the committed corpus or
        tests/*.m (fast; used to verify the emitters before the full run).
    sage -python tests/generate_corpus.py --inputs-only   # re-runs assembly
        only and writes tests/corpus_curves.json; no oracle computation, no
        `expect` blocks, no tests/*.m (fast; use after editing assemble_inputs
        to refresh the corpus ahead of the full run)."""
import datetime
import json
import os
import sys
from sage.all import (QQ, ZZ, PolynomialRing, NumberField, EllipticCurve,
                      kronecker_symbol, set_random_seed, prime_range, gcd)
from sage.arith.misc import fundamental_discriminant
from sage.parallel.decorate import fork, parallel
from sage.parallel.ncpus import ncpus as _sage_ncpus
try:
    from sage.schemes.elliptic_curves.mod_poly import classical_modular_polynomial
except Exception:      # optional; the modular-polynomial screen degrades to skip
    classical_modular_polynomial = None

SEED = 20260720
MAZUR = [2, 3, 5, 7, 11, 13, 17, 19, 37, 43, 67, 163]
ORACLE_TIMEOUT = 120   # cap (seconds) per oracle, and per-ell for the oE
                       # construction certificate (steps 2-3 of the cascade)
ENTRY_TIMEOUT = 1800   # outer per-entry fork cap (compute_oracles): a backstop
                       # above ORACLE_TIMEOUT in case several per-oracle/per-ell
                       # caps chain up past it on one curve; not exercised by the
                       # current corpus (see the "entry_timeout" drop marker)
FROB_SCREEN_PRIMES = 40   # good primes cached per curve for the Frobenius screen
_ACTIVE_CAP = ORACLE_TIMEOUT   # set by compute_oracles(); read by capped()

def _worker_count():
    """Entry-level fork width for compute_oracles: min(cores - 2, 16) floored at
    4. sage.parallel.ncpus.ncpus() is the spec-named source, but it honors
    SAGE_NUM_THREADS, which this sandbox pins to 1 (a per-process BLAS thread cap,
    set alongside OPENBLAS_NUM_THREADS=1) -- that collapses ncpus() to 1, which is
    a library-thread budget, not a fork budget. Fall back to the CPU affinity when
    ncpus() is degenerate so the entry parallelism actually spreads across cores."""
    n = _sage_ncpus()
    if n <= 1:
        aff = os.sched_getaffinity(0) if hasattr(os, "sched_getaffinity") else None
        n = len(aff) if aff else (os.cpu_count() or 1)
    return max(4, min(n - 2, 16))

NCPUS = _worker_count()

def K_from(coeffs):
    R = PolynomialRing(QQ, 'x')
    if coeffs == [0, 1]:
        return QQ
    return NumberField(R(coeffs), 'w')

def vec(K, elt):
    return [QQ(c) for c in (list(K(elt)) if K is not QQ else [elt])]

def entry(id_, stratum, K, E, E2=None):
    e = {"id": id_, "stratum": stratum,
         "field": [QQ(c) for c in (K.defining_polynomial().list() if K is not QQ else [0,1])],
         "curve": {"model": "ainvs", "ainvs": [vec(K, a) for a in E.ainvs()]}}
    if E2 is not None:
        e["curve2"] = {"model": "ainvs", "ainvs": [vec(K, a) for a in E2.ainvs()]}
    return e

def assemble_inputs(data):
    # Idempotent: drop every previously-assembled entry before appending fresh
    # ones, so re-running assembly regenerates exactly the fixture/differential
    # strata instead of duplicating them. The 250 lmfdb-qcurve entries (Task 1)
    # and _provenance are untouched.
    data["entries"] = [e for e in data["entries"] if e.get("stratum") == "lmfdb-qcurve"]
    es = data["entries"]
    set_random_seed(SEED)
    R = PolynomialRing(QQ, 'x'); x = R.gen()
    # fixture-ex58: Billerey Example 5.8, K = Q(sqrt-3, sqrt-7), conductor 2*O_K.
    K58 = NumberField((x**2 + 3, x**2 + 7), names=('s3', 's7')).absolute_field('t')
    r3 = (x**2 + 3).roots(K58, multiplicities=False)[0]
    r7 = (x**2 + 7).roots(K58, multiplicities=False)[0]
    r21 = r3 * r7          # sqrt(-3)*sqrt(-7) = sqrt(21)
    a4 = QQ(81)/4 * (69 + 43*r3 + 29*r7 + 17*r21)
    a6 = 162 * (207 - 84*r3 - 54*r7 + 46*r21)
    E58 = EllipticCurve(K58, [0, 0, 0, a4, a6])
    assert not E58.has_cm()
    es.append(entry("fixture-ex58", "fixture-ex58", K58, E58))
    # gate fixtures over Q(sqrt 29): 11a1 base-changed; 3 is inert, choose a
    # split prime too (7 splits since kronecker(29,7)=1).
    K29 = NumberField(x**2 - 29, 'w')
    E11 = EllipticCurve(QQ, [0, -1, 1, -10, -20]).change_ring(K29)
    es.append(entry("fixture-gate-inert", "fixture-gate-inert", K29, E11))
    es.append(entry("fixture-gate-split", "fixture-gate-split", K29, E11))
    # fixture-ordq: field with class number > 1 and a prime with ord([q]) < h;
    # K = Q(sqrt -23) has h = 3; also record a nonprincipal prime exists.
    K23 = NumberField(x**2 + 23, 'w')
    es.append(entry("fixture-ordq", "fixture-ordq", K23,
                    EllipticCurve(QQ, [1, 1, 0, 1, 0]).change_ring(K23)))
    # fixture-deg1-fldnum: degree-one FldNum must dispatch to branch 1.
    K1 = NumberField(x - 1, 'w')
    es.append(entry("fixture-deg1-fldnum", "fixture-deg1-fldnum", K1,
                    EllipticCurve(K1, [0, -1, 1, -10, -20])))
    # fixture-localglobal: j = 2268945/128 over Q(i) (7 passes local tests, no
    # global 7-isogeny; verified in review round 3).
    Ki = NumberField(x**2 + 1, 'i')
    Elg = EllipticCurve(Ki, j=QQ(2268945)/128)
    es.append(entry("fixture-localglobal", "fixture-localglobal", Ki, Elg))
    # CM fixtures (spec fixture 6):
    es.append(entry("fixture-cm-1728-Qi", "fixture-cm", Ki, EllipticCurve(Ki, [1, 0])))
    es.append(entry("fixture-cm-0-Qi", "fixture-cm", Ki, EllipticCurve(Ki, [0, 1])))
    es.append(entry("fixture-cm-1728-Q", "fixture-cm", QQ, EllipticCurve(QQ, [1, 0])))
    K5 = NumberField(x**2 - 5, 'w')
    es.append(entry("fixture-cm-notinK", "fixture-cm", K5,
                    EllipticCurve(QQ, [0, 1]).change_ring(K5)))   # CM by Q(sqrt-3), not in K
    # non-maximal order: j = 66^3 has CM by disc -16 = f^2 * (-4), f = 2.
    es.append(entry("fixture-cm-nonmax", "fixture-cm", Ki, EllipticCurve(Ki, j=QQ(66)**3)))
    # congruence fixtures (spec fixture 8):
    E14, E14b = EllipticCurve(QQ, '14a1'), EllipticCurve(QQ, '14a2')
    es.append(entry("fixture-cong-identical", "fixture-cong-isogenous", QQ, E14, E14))
    es.append(entry("fixture-cong-isogpair", "fixture-cong-isogenous", QQ, E14, E14b))
    Etw = E14.quadratic_twist(5)
    es.append(entry("fixture-cong-twist", "fixture-cong-twist", QQ, E14, Etw))
    # ell=2 congruence: 15a1 and 17a2 both have full rational 2-torsion
    # (two_torsion_rank 2), so Galois acts trivially on E[2] for both, giving
    # a_q(E1) - a_q(E2) = 0 mod 4 for every good q (injectivity of torsion
    # under reduction); different conductors rule out isogeny. Verified: the
    # trace-gcd oracle gives primes = [2] (gcd = 4 over q < 8000, stabilizing
    # at bound 2000). The former pair (11a1, 14a1) was NOT 2-congruent (oracle
    # gcd = 1, i.e. empty primes) despite the fixture's name.
    es.append(entry("fixture-cong-2", "fixture-cong-finite", QQ,
                    EllipticCurve(QQ, '15a1'), EllipticCurve(QQ, '17a2')))
    es.append(entry("fixture-cong-control", "fixture-cong-finite", QQ,
                    EllipticCurve(QQ, '11a1'), EllipticCurve(QQ, '37a1')))
    # sage-doc examples over number fields (known reducible_primes docstrings):
    K5r = NumberField(x**2 - 29, 'a')
    es.append(entry("sage-doc-29", "sage-doc", K5r,
                    EllipticCurve(K5r, [1, 0, ((5 + K5r.gen())/2)**2, 0, 0])))
    # base changes (spec fixture 10):
    for lab in ['11a1', '37a1', '389a1']:
        E = EllipticCurve(QQ, lab)
        es.append(entry(f"basechange-{lab}-Q", "basechange-q", QQ, E))
        es.append(entry(f"basechange-{lab}-Qs29", "basechange-K", K29, E.change_ring(K29)))
    # differential strata, seeded:
    for i in range(10):   # generic controls
        for Kd, tag in [(K29, 'q29'), (Ki, 'qi')]:
            a = [Kd([ZZ.random_element(-9, 10) for _ in range(2)]) for _ in range(2)]
            try:
                E = EllipticCurve(Kd, [a[0], a[1]])
            except ArithmeticError:
                continue
            es.append(entry(f"diff-generic-{tag}-{i}", "diff-generic", Kd, E))
    # X_0(ell) stratum: curves over K with a genuine ell-isogeny, via the
    # classical rational parametrization of the j-line by X_0(ell), ell in
    # {2,3,5,7,13} (genus 0): pick random t in K, build E with j = j_ell(t).
    # j_2 and j_13 below are Sage's own Fricke_module(ell) (in
    # isogeny_small_degree.py): the brief's original j_2 = 256(x^2+x+16)^3 /
    # (x^2(x+16)) and j_13 = .../x^13 do not land on X_0(2)/X_0(13) at all
    # (isogenies_prime_degree came back empty for essentially every probed t,
    # over both Q and K29); j_3, j_5, j_7 passed the same probe and are kept.
    J = {2: (x+16)**3 / x,
         3: ((x+27)*(x+3)**3)/x,
         5: ((x**2+250*x+3125)**3)/x**5,
         7: ((x**2+13*x+49)*(x**2+245*x+2401)**3)/x**7,
         13: ((x**2+5*x+13)*(x**4+7*x**3+20*x**2+19*x+1)**3)/x}
    for ell, jmap in J.items():
        for i in range(3):
            t = K29([ZZ.random_element(1, 50), ZZ.random_element(0, 3)])
            jv = jmap(t)
            if jv in [0, 1728]:
                continue
            E = EllipticCurve(K29, j=jv)
            # verify the parametrization actually delivers an ell-isogeny
            # before trusting this as a diff-x0 (and possible diff-congpair)
            # ground-truth entry.
            phis = E.isogenies_prime_degree(ell)
            if not phis:
                continue
            es.append(entry(f"diff-x0-{ell}-{i}", "diff-x0", K29, E))
            if i == 0:   # congruence pair stratum: E and its ell-isogenous partner
                es.append(entry(f"diff-congpair-{ell}", "diff-congpair",
                                K29, E, phis[0].codomain()))
    # CM with j not in Q: disc -15 (h=2), j in Q(sqrt 5).
    from sage.schemes.elliptic_curves.cm import hilbert_class_polynomial
    H = hilbert_class_polynomial(-15)
    jr = H.roots(K5, multiplicities=False)
    for i, jv in enumerate(jr):
        es.append(entry(f"diff-cmj-15-{i}", "diff-cmj", K5, EllipticCurve(K5, j=jv)))
    return data

# ===========================================================================
# Task 3: oracle computation and Magma test-file emitters.
#
# The oracles below reconstruct each entry's K and E in Sage the same way the
# emitted Magma BuildField/BuildCurve preamble does. Sage and Magma may pick
# conjugate roots of the defining polynomial, but every oracle value used in an
# assert (O(E), CM data, congruence-trace gcd) is a Galois invariant, and the
# golden charpolys are recorded for a Q-rational curve, so the reconstruction is
# consistent on both sides.
# ===========================================================================

def sage_field(coeffs):
    R = PolynomialRing(QQ, 'x')
    cc = [QQ(c) for c in coeffs]
    if cc == [QQ(0), QQ(1)]:
        return QQ
    return NumberField(R(cc), 'w')

def sage_elt(K, vecc):
    v = [QQ(c) for c in vecc]
    return v[0] if K is QQ else K(v)

def sage_curve(K, model):
    if model["model"] == "AB":
        return EllipticCurve(K, [0, 0, 0, sage_elt(K, model["A"]),
                                 sage_elt(K, model["B"])])
    return EllipticCurve(K, [sage_elt(K, a) for a in model["ainvs"]])

def rebuild(e):
    K = sage_field(e["field"])
    E = sage_curve(K, e["curve"])
    E2 = sage_curve(K, e["curve2"]) if "curve2" in e else None
    return K, E, E2

def is_deg1(K):
    return (K is QQ) or (K.absolute_degree() == 1)

def capped(name, thunk, dropped):
    """Run thunk in a forked subprocess with a hard timeout (_ACTIVE_CAP) and
    return its result. fork() kills the subprocess on timeout, so this bounds even
    the PARI/C sections of isogenies_prime_degree that alarm()/AlarmInterrupt
    cannot interrupt (a huge-coefficient X_0 curve hangs an alarm for minutes).
    On timeout, or an unhandled error in the child, fork returns a 'NO DATA'
    sentinel string; record the drop and return None."""
    runner = fork(lambda _: thunk(), timeout=_ACTIVE_CAP, verbose=False)
    r = runner(0)
    if isinstance(r, str) and r.startswith("NO DATA"):
        dropped.append(name)
        return None
    return r

def deg1_curve_over_Q(e):
    """The isomorphic Q-model of a degree-one entry (a degree-1 NumberField is
    canonically Q). Used so O(E) is computed with fast rational arithmetic:
    isogenies_prime_degree over a degree-1 NumberField times out."""
    m = e["curve"]
    if m["model"] == "AB":
        return EllipticCurve(QQ, [0, 0, 0, QQ(m["A"][0]), QQ(m["B"][0])])
    return EllipticCurve(QQ, [QQ(a[0]) for a in m["ainvs"]])

# --- oE screens (cheapest-certain-first cascade) ----------------------------
#
# DIRECTION DISCIPLINE: the two screens below only ever produce CERTAIN
# NEGATIVES (ell provably NOT in O(E)); membership in O(E) is established ONLY by
# the construction certificate (step 3). So the screened O(E) equals the
# brute-force isogenies_prime_degree sweep exactly -- a screen can only remove an
# ell that construction would also have found absent -- while running orders of
# magnitude fewer of the expensive constructor calls.

def _quad_irreducible_mod_ell(a, n, ell):
    """Is x^2 - a x + n irreducible over GF(ell)? (a = trace, n = norm of a good
    Frobenius.) ell = 2: only x^2+x+1 is irreducible, i.e. a and n both odd.
    ell odd: irreducible iff the discriminant a^2 - 4n is a nonzero non-square."""
    if ell == 2:
        return (n % 2 == 1) and (a % 2 == 1)
    D = (a * a - 4 * n) % ell
    if D == 0:
        return False
    return kronecker_symbol(D, ell) == -1

def _frob_data(E, base, n_primes=FROB_SCREEN_PRIMES):
    """Cache (residue_char, a_q, Norm(q)) for the n_primes good primes of E of
    smallest norm, reused by the Frobenius screen across every ell. E.reduction
    raises at bad primes, so try/except keeps only good ones -- no conductor
    computation needed. Over QQ the primes are p; over K they are prime ideals
    taken norm-ascending. Runs under the caller's capped() (up to n_primes
    PARI-backed E.reduction(...).trace_of_frobenius() calls, otherwise
    unbounded)."""
    data = []
    if base is QQ:
        for p in prime_range(2, 2000):
            try:
                a = ZZ(E.reduction(p).trace_of_frobenius())
            except Exception:
                continue
            data.append((int(p), int(a), int(p)))
            if len(data) >= n_primes:
                break
    else:
        cand = [(int(P.norm()), int(P.smallest_integer()), P)
                for P in base.primes_of_bounded_norm(500)]
        cand.sort(key=lambda t: t[0])
        for nrm, p, P in cand:
            try:
                a = ZZ(E.reduction(P).trace_of_frobenius())
            except Exception:
                continue
            data.append((int(p), int(a), int(nrm)))
            if len(data) >= n_primes:
                break
    return data

def _frob_screen_kills(frob, ell):
    """CERTAIN NEGATIVE: True if some cached good prime q (residue char != ell)
    has Frobenius charpoly irreducible mod ell, proving E[ell] is irreducible
    over K (a K-rational ell-subgroup would give Frobenius an eigenvalue mod ell,
    hence a charpoly root, at every good q coprime to ell)."""
    for (res_char, aq, nq) in frob:
        if res_char == ell:
            continue
        if _quad_irreducible_mod_ell(aq, nq, ell):
            return True
    return False

def _modpoly_no_root(base, jE, ell):
    """CERTAIN NEGATIVE: True if Phi_ell(j(E), Y) has no root in the base field,
    so E has no ell-isogeny over K (a K-rational ell-subgroup C gives the
    K-rational root j(E/C)). Caller guards j not in {0, 1728} and swallows errors
    (the screen then just skips, never a false drop)."""
    Phi = classical_modular_polynomial(int(ell))
    d = Phi.dict()
    max_ex = max(ex for (ex, ey) in d)
    jb = base(jE)
    jpow = [base(1)] * (max_ex + 1)
    for k in range(1, max_ex + 1):
        jpow[k] = jpow[k - 1] * jb
    coeffs = {}
    for (ex, ey), c in d.items():
        coeffs[ey] = coeffs.get(ey, base(0)) + base(c) * jpow[ex]
    f = PolynomialRing(base, 'Y')(coeffs)
    return not f.roots(base, multiplicities=False)

def oracle_oE(E, base, deg1, dropped, stats):
    """O(E) = { ell : isogenies_prime_degree(ell) nonempty }, decided by the
    cheapest-certain-first cascade (Frobenius screen, then modular-polynomial
    screen, then the construction certificate). Full Mazur list (incl. 163) over
    degree-one fields; primes <= 100 over number fields. E is the curve to sweep
    with base ring `base` (a Q-model for degree-one entries; see
    deg1_curve_over_Q). _frob_data itself runs under the ORACLE_TIMEOUT fork cap;
    a timeout records an "oE:frob" drop and degrades to an empty screen (never a
    false certain-negative -- see DIRECTION DISCIPLINE above: an unavailable
    screen just contributes nothing, falling through to modpoly/construction for
    every ell). Each surviving ell's steps 2-3 run under their own per-ell
    ORACLE_TIMEOUT fork cap; a timeout records the (curve, ell) drop, never
    silently. `stats` is filled with the per-screen hit counts for logging."""
    sweep = [ZZ(p) for p in MAZUR] if deg1 else [ZZ(p) for p in prime_range(101)]
    frob = capped("oE:frob", lambda: _frob_data(E, base), dropped)
    if frob is None:   # timed out: no screen data, never a false negative
        frob = []
    jE = E.j_invariant()
    j_special = jE in (base(0), base(1728))   # skip modpoly screen at 0/1728
    out = []
    n_frob = n_modpoly = n_in = n_out = n_drop = 0
    for ell in sweep:
        elli = int(ell)
        # Step 1: Frobenius screen -- fast, inline, certain negative.
        if _frob_screen_kills(frob, elli):
            n_frob += 1
            continue
        do_modpoly = (not j_special) and (classical_modular_polynomial is not None)
        def per_ell(ell=ell, do_modpoly=do_modpoly):
            # ell must stay a Sage Integer: a Python int trips
            # isogenies_prime_degree's internal l.is_prime() call.
            try:
                # Step 2: modular-polynomial screen -- certain negative.
                if do_modpoly:
                    try:
                        if _modpoly_no_root(base, jE, ell):
                            return "modpoly_drop"
                    except Exception:
                        pass   # screen inconclusive; fall through to construction
                # Step 3: construction certificate -- the only positive path.
                # Over number fields pass minimal_models=False (we need only
                # nonemptiness; building codomain minimal models is wasted work
                # and trips a Sage round() bug on some inputs). Over QQ that kwarg
                # does not exist.
                if deg1:
                    phis = E.isogenies_prime_degree(ell)
                else:
                    phis = E.isogenies_prime_degree(ell, minimal_models=False)
                return "in" if phis else "out"
            except Exception:
                # No constructible isogeny or a Sage arithmetic failure for this
                # ell: certainly not certified in. Returning here (rather than
                # raising) keeps the fork's 'NO DATA' sentinel meaning ONLY a
                # hard timeout/crash, i.e. the drop case.
                return "out"
        r = capped("oE:ell=%d" % elli, per_ell, dropped)
        if r is None:                 # timeout: capped() already recorded the drop
            n_drop += 1
        elif r == "in":
            n_in += 1
            out.append(elli)
        elif r == "modpoly_drop":
            n_modpoly += 1
        else:
            n_out += 1
    stats.update(sweep=len(sweep), frob=n_frob, modpoly=n_modpoly,
                 construct_in=n_in, construct_out=n_out, dropped=n_drop)
    return out

def oracle_cm(E, K, dropped):
    """Geometric CM (E.has_cm() over Qbar). cm_discriminant() is the ORDER
    discriminant D_O = f^2 D_F; derive D_F and f, and CMInBaseField as the
    IsSquare test on D_F (never square over Q, so degree-1 CM is not in K)."""
    def go():
        if not E.has_cm():
            return {"is_cm": False}
        DO = ZZ(E.cm_discriminant())
        DF = ZZ(fundamental_discriminant(DO))
        f = ZZ(DO // DF).isqrt()
        ib = bool(K(DF).is_square()) if K is not QQ else bool(QQ(DF).is_square())
        return {"is_cm": True, "order_disc": int(DO), "fund_disc": int(DF),
                "f": int(f), "in_base_field": ib}
    return capped("cm", go, dropped)

def oracle_soft(E, dropped):
    """Sage reducible_primes (SOFT, print-only tightness comparison)."""
    def go():
        try:
            rp = E.galois_representation().reducible_primes()
        except Exception:   # print-only; AlarmInterrupt (drop) still propagates
            return None
        return sorted(int(p) for p in rp)
    return capped("soft", go, dropped)

def oracle_golden(E, K, dropped):
    """Charpoly triples for the gate curve over Q(sqrt29): p=7 splits, p=3 is
    inert, both good for 11a1. E is Q-rational, so the trace is prime-above
    independent and the recorded [Norm, -a, 1] matches any prime above p."""
    def go():
        res = {}
        for tag, p in (("split", 7), ("inert", 3)):
            P = K.primes_above(p)[0]
            fp = E.reduction(P).frobenius_polynomial()
            res[tag] = {"p": int(p), "norm": int(P.norm()),
                        "charpoly": [int(c) for c in fp.list()]}
        return res
    return capped("golden", go, dropped)

def _cong_terms(E1, E2, K, max_norm):
    """[(Norm(q), p_q*(a_q(E1)-a_q(E2)))] over prime ideals q good for both
    curves (coprime to both conductors), Norm(q) <= max_norm. The gcd of the
    second components is order-independent, so no sorting is needed here."""
    terms = []
    if K is QQ:
        N1, N2 = ZZ(E1.conductor()), ZZ(E2.conductor())
        for p in prime_range(max_norm + 1):
            if N1 % p == 0 or N2 % p == 0:
                continue
            terms.append((int(p), ZZ(p) * (ZZ(E1.ap(p)) - ZZ(E2.ap(p)))))
    else:
        N1, N2 = E1.conductor(), E2.conductor()
        for P in K.primes_of_bounded_norm(max_norm):
            if P.divides(N1) or P.divides(N2):
                continue
            p = ZZ(P.smallest_integer())   # rational prime below P (residue char)
            a1 = ZZ(E1.reduction(P).trace_of_frobenius())
            a2 = ZZ(E2.reduction(P).trace_of_frobenius())
            terms.append((int(P.norm()), p * (a1 - a2)))
    return terms

def _cong_escalate(terms, norm_bound, max_norm_bound):
    """Spec-pinned B-phase escalation on G := gcd of the terms: zero-sentinel
    (G == 0 is a formal TOP, never a plateau, never factored), doubling from
    norm_bound to max_norm_bound, plateau stop (a set unchanged across a doubling,
    including one landing exactly at the cap), G == 1 absorbing early exit, and
    the single-evaluation rule (stabilized iff the set is empty). Returns the
    reproduction expect {primes, kind, stabilized, bound}: kind 'Finite' if the
    final G is nonzero, else 'ZeroAtCap' (the engine then goes to the BFS at
    degree >= 2 or returns Undecided at degree 1)."""
    def gcd_upto(b):
        G = ZZ(0)
        for (nrm, t) in terms:
            if nrm <= b:
                G = gcd(G, t)
        return G
    bound = norm_bound
    prev = None
    while True:
        G = gcd_upto(bound)
        if G == 1:
            return {"primes": [], "kind": "Finite", "stabilized": True,
                    "bound": int(bound)}
        cur = None if G == 0 else frozenset(int(p) for p in G.prime_divisors())
        if prev is not None and cur is not None and cur == prev:
            return {"primes": sorted(cur), "kind": "Finite",
                    "stabilized": True, "bound": int(bound)}
        if bound >= max_norm_bound:
            if cur is None:
                return {"primes": None, "kind": "ZeroAtCap",
                        "stabilized": False, "bound": int(bound)}
            return {"primes": sorted(cur), "kind": "Finite",
                    "stabilized": False, "bound": int(bound)}
        prev = cur
        bound = min(2 * bound, max_norm_bound)

def oracle_cong(E1, E2, K, dropped, norm_bound=1000, max_norm_bound=8000,
                name="cong"):
    def go():
        terms = _cong_terms(E1, E2, K, max_norm_bound)
        return _cong_escalate(terms, norm_bound, max_norm_bound)
    return capped(name, go, dropped)

GATE_IDS = {"fixture-gate-inert", "fixture-gate-split"}

def _compute_expect_for_entry(e):
    """Pure function of one entry: rebuild K/E and return (id, expect, oE_stats).
    Runs inside an entry-level fork; the per-oracle capped() calls fork again
    (nested), and every returned value is a plain JSON type so the parent can
    pickle it back and reattach by id. `expect["dropped"]` lists this entry's
    capped drops (oracle name, or "oE:ell=<ell>" for a per-ell construction)."""
    eid = e["id"]
    dropped = []
    oe_stats = {}
    K, E, E2 = rebuild(e)
    expect = {"deg1": is_deg1(K)}
    if E2 is not None:
        # Only the fixture-cong-* pairs are asserted (congruence section);
        # diff-congpair pairs are differential-only (the driver calls the
        # engine and prints; no oracle) and are isogenous, so the engine
        # returns AllPrimes, not the gcd's Finite.
        if e["stratum"].startswith("fixture-cong"):
            expect["cong"] = oracle_cong(E, E2, K, dropped)
            if e["stratum"] == "fixture-cong-twist":
                expect["cong_lowbound"] = oracle_cong(
                    E, E2, K, dropped, norm_bound=2, max_norm_bound=2,
                    name="cong_lowbound")
    else:
        cm = oracle_cm(E, K, dropped)
        expect["cm"] = cm
        is_cm = bool(cm and cm.get("is_cm"))
        Eo = deg1_curve_over_Q(e) if (expect["deg1"] and K is not QQ) else E
        expect["oE"] = oracle_oE(Eo, Eo.base_ring(), expect["deg1"], dropped,
                                 oe_stats)
        if (not expect["deg1"]) and (not is_cm):
            expect["soft_reducible"] = oracle_soft(E, dropped)
        if eid in GATE_IDS:
            expect["golden"] = oracle_golden(E, K, dropped)
    expect["dropped"] = dropped
    return (eid, expect, oe_stats)

def compute_oracles(entries, cap=ORACLE_TIMEOUT, entry_timeout=ENTRY_TIMEOUT):
    """Attach an `expect` block to every entry in place, computing entries in
    parallel over NCPUS forks, each capped at `entry_timeout` wall-clock seconds
    -- a backstop above the per-oracle `cap` (e.g. many per-ell oE drops chaining
    up on one curve could in principle exceed the per-oracle cap's sum without
    it). Returns the aggregated [(id, oracle)] drop list for the emitted header.
    Each entry's expect block is a pure function of the entry, so results are
    reattached BY ID in the original entry order: the written corpus and emitted
    .m files are byte-identical regardless of the order the forks finish. A
    whole-entry timeout (or crash) is recorded as an "entry_timeout" drop rather
    than crashing compute_oracles or silently omitting the entry: the entry
    keeps its id and gets a minimal expect block, so the emitters' existing
    None-handling (oE/cm/etc. all absent) skips its asserts cleanly instead of
    asserting something wrong."""
    global _ACTIVE_CAP
    _ACTIVE_CAP = cap
    n = len(entries)
    worker = parallel(ncpus=NCPUS, timeout=entry_timeout)(_compute_expect_for_entry)
    results = {}
    tot = {}
    done = 0
    for (_args, out) in worker([(e,) for e in entries]):
        if isinstance(out, str) and out.startswith("NO DATA"):
            # Whole-entry fork timeout/crash: recover the entry from the
            # parent-side input (never round-tripped through the dead child,
            # so it's available even though nothing came back). deg1 is cheap
            # and safe to recompute here -- pure NumberField-degree check on
            # stored coefficients, no PARI reduction/isogeny calls -- so the
            # emitters' branch1-vs-branch2 classification still matches the
            # entry's real shape; every other field is left absent, which the
            # emitters already treat as "oracle dropped" (see e.g. `oE is
            # None` in _sec_branch1/_sec_branch2/_sec_cm).
            e = _args[0][0]
            eid = e["id"]
            try:
                deg1 = is_deg1(sage_field(e["field"]))
            except Exception:
                deg1 = None
            expect = {"deg1": deg1, "dropped": ["entry_timeout"]}
        else:
            eid, expect, oe_stats = out
            for k, v in oe_stats.items():
                tot[k] = tot.get(k, 0) + v
        results[eid] = expect
        done += 1
        drp = expect.get("dropped") or []
        print("oracle %d/%d %s%s" % (done, n, eid,
              ("  DROPPED:" + ",".join(drp)) if drp else ""), flush=True)
    # Reattach in the original entry order so the output is deterministic.
    agg = []
    for e in entries:
        expect = results[e["id"]]
        e["expect"] = expect
        agg.extend((e["id"], nm) for nm in (expect.get("dropped") or []))
    if tot:
        print("oE screens (totals over swept entries): sweep=%d frob_killed=%d "
              "modpoly_killed=%d construct_in=%d construct_out=%d dropped=%d"
              % (tot.get("sweep", 0), tot.get("frob", 0), tot.get("modpoly", 0),
                 tot.get("construct_in", 0), tot.get("construct_out", 0),
                 tot.get("dropped", 0)), flush=True)
    return agg

# --- Magma literal helpers ---------------------------------------------------

def _mnum(v):
    return str(v)   # JSON ints/strings are valid Magma integer/rational literals

def _mseq(coords):
    return "[" + ", ".join(_mnum(c) for c in coords) + "]"

def _ainvs_rows(model):
    if model["model"] == "AB":
        return [[0], [0], [0], model["A"], model["B"]]
    return model["ainvs"]

def _mainvs(model):
    return "[ " + ", ".join(_mseq(r) for r in _ainvs_rows(model)) + " ]"

def _mfield(coeffs):
    return "[" + ", ".join(_mnum(c) for c in coeffs) + "]"

def _mset(primes):
    return ("{ " + ", ".join(str(p) for p in primes) + " }") if primes else "{ Integers() | }"

def _mlist(ints):
    return ("[ " + ", ".join(str(p) for p in ints) + " ]") if ints else "[ Integers() | ]"

def _build1(e):
    return ["    K := BuildField(%s);" % _mfield(e["field"]),
            "    E := BuildCurve(K, %s);" % _mainvs(e["curve"])]

def _build2(e):
    return ["    K := BuildField(%s);" % _mfield(e["field"]),
            "    E1 := BuildCurve(K, %s);" % _mainvs(e["curve"]),
            "    E2 := BuildCurve(K, %s);" % _mainvs(e["curve2"])]

PREAMBLE = r'''// ==========================================================================
// GENERATED by tests/generate_corpus.py -- DO NOT EDIT BY HAND.
%(header)s// ==========================================================================
if assigned spec then useSpec := spec; else useSpec := "magma/spec"; end if;
if useSpec ne "" then
    // CHIMP supplies the star/bracket charpoly intrinsics the engine needs.
    ok, _ := IsIntrinsic("PowerCharacteristicPolynomial");
    error if not ok, "CHIMP is not attached; AttachSpec your CHIMP.spec first";
    AttachSpec(useSpec);
end if;
// useSpec eq "" is the red-state syntactic check: no engine spec, CHIMP not needed.
if not assigned section then section := "all"; end if;
if not assigned cmscope then cmscope := "1"; end if;

R<x> := PolynomialRing(Rationals());

function BuildField(coeffs)
    if coeffs eq [0, 1] then return Rationals(); end if;
    // NumberField(R!coeffs) collapses a degree-1 defining polynomial to
    // FldRat in this Magma version; force a genuine FldNum here so the
    // degree-one dispatch path is actually exercised.
    if #coeffs eq 2 then return RationalsAsNumberField(); end if;
    return NumberField(R ! coeffs);
end function;

function Elt(K, v)
    if Type(K) eq FldRat then return K ! v[1]; end if;
    return &+[ (Rationals() ! v[i]) * K.1^(i-1) : i in [1..#v] ];
end function;

function BuildCurve(K, ainvs)
    return EllipticCurve([ Elt(K, a) : a in ainvs ]);
end function;

function PrimeAbove(K, p, i)
    return Decomposition(Integers(K), p)[i][1];
end function;
'''

def _sec_golden(gates):
    L = ["procedure Test_golden()"]
    for e in gates:
        eid = e["id"]
        g = e["expect"].get("golden")
        if g is None:
            L.append('    printf "  SKIP %o: golden dropped\\n", "' + eid + '";')
            continue
        L += _build1(e)
        sp, ip = g["split"], g["inert"]
        L.append("    qs := PrimeAbove(K, %d, 1);" % sp["p"])
        L.append("    qi := PrimeAbove(K, %d, 1);" % ip["p"])
        L.append("    assert [Integers() ! c : c in Coefficients(FrobeniusCharpoly(E, qs))] eq %s;  // %s split charpoly" % (_mlist(sp["charpoly"]), eid))
        L.append("    assert [Integers() ! c : c in Coefficients(FrobeniusCharpoly(E, qi))] eq %s;  // %s inert charpoly" % (_mlist(ip["charpoly"]), eid))
        if eid == "fixture-gate-inert":
            L.append("    assert BillereyRq(E, qi) eq BillereyBl(E, %d);  // inert gate R_q = B_l" % ip["p"])
        else:
            L.append("    assert (BillereyRq(E, qs) ne BillereyBl(E, %d)) or (BillereyBl(E, %d) eq 0);  // split: gate must not hold" % (sp["p"], sp["p"]))
    L.append('    printf "SECTION golden: PASS\\n";')
    L.append("end procedure;\n")
    return "\n".join(L)

def _sec_branch1(entries):
    L = ["procedure Test_branch1()"]
    for e in entries:
        eid = e["id"]; oE = e["expect"].get("oE")
        L += _build1(e)
        L.append("    L, info := IsogenyPrimes(E);")
        L.append('    assert info`Source eq "IsogenyPrimes";  // %s' % eid)
        L.append('    assert info`Kind eq "Finite";')
        L.append("    assert info`Exact;")
        if oE is None:
            L.append('    printf "  SKIP %o: oE dropped\\n", "' + eid + '";')
        else:
            L.append("    assert Set(L) eq %s;  // exact O(E), complete by Mazur" % _mset(oE))
            L.append("    for ell in %s do assert MayBeReducible(ell, L, info); end for;" % _mlist(oE))
    L.append('    printf "SECTION branch1: PASS\\n";')
    L.append("end procedure;\n")
    return "\n".join(L)

def _sec_branch2(entries):
    L = ["procedure Test_branch2()"]
    for e in entries:
        eid = e["id"]; exp = e["expect"]; oE = exp.get("oE"); soft = exp.get("soft_reducible")
        L += _build1(e)
        L.append("    L, info := IsogenyPrimes(E);")
        L.append('    assert info`Source eq "IsogenyPrimes";  // %s' % eid)
        L.append("    assert not info`Exact;")
        if oE is None:
            L.append('    printf "  SKIP %o: oE dropped\\n", "' + eid + '";')
        else:
            L.append("    for ell in %s do assert MayBeReducible(ell, L, info); end for;  // containment" % _mlist(oE))
        if eid == "fixture-localglobal":
            L.append("    assert 7 in Set(L);            // local-global false positive kept")
            if oE is not None:
                L.append("    assert 7 notin %s;            // ... but no global 7-isogeny" % _mset(oE))
        if eid == "fixture-ex58":
            L.append("    assert info`BoundsUsed[2] ne 0;  // R-phase ran (B_l = 0 for all l)")
        if soft is None:
            L.append('    printf "  branch2 %o soft=UNAVAILABLE L=%o\\n", "' + eid + '", Sort(SetToSequence(Set(L)));')
        else:
            L.append('    printf "  branch2 %o soft=%o L=%o\\n", "' + eid + '", ' + _mlist(soft) + ', Sort(SetToSequence(Set(L)));')
    L.append('    printf "SECTION branch2: PASS\\n";')
    L.append("end procedure;\n")
    return "\n".join(L)

def _sec_cm(entries):
    L = ["procedure Test_cm()"]
    for e in entries:
        eid = e["id"]; exp = e["expect"]; oE = exp.get("oE"); cm = exp.get("cm") or {}
        in_bf = bool(cm.get("in_base_field"))
        kind = "CMFamily" if in_bf else "Finite"
        exact = bool(exp.get("deg1"))   # degree-1 CM takes branch 1 (Exact true)
        L += _build1(e)
        L.append("    L, info := IsogenyPrimes(E);")
        L.append('    assert info`Source eq "IsogenyPrimes";  // %s' % eid)
        L.append("    assert info`IsCM;")
        L.append("    assert info`CMOrderDiscriminant eq %d;" % cm["order_disc"])
        L.append("    assert info`CMFundamentalDiscriminant eq %d;" % cm["fund_disc"])
        L.append("    assert info`CMConductor eq %d;" % cm["f"])
        L.append("    assert info`CMInBaseField eq %s;" % ("true" if in_bf else "false"))
        L.append('    assert info`Kind eq "%s";' % kind)
        L.append("    assert %sinfo`Exact;" % ("" if exact else "not "))
        if oE is None:
            L.append('    printf "  SKIP %o: oE dropped\\n", "' + eid + '";')
        elif kind == "Finite":
            L.append("    assert %s subset Set(L);  // O(E) subset L" % _mset(oE))
        else:
            L.append("    for ell in %s do assert MayBeReducible(ell, L, info); end for;  // CMFamily denotation" % _mlist(oE))
        if eid == "fixture-cm-1728-Qi":
            if oE is not None:
                L.append("    assert 13 in %s;  // 13 construction-reducible (F-family)" % _mset(oE))
            L.append("    assert HasPrimeIsogeny(E, 5);  // split sample: F subset R(E)")
            L.append("    assert HasPrimeIsogeny(E, 2);  // ramified sample")
        if eid == "fixture-cm-nonmax":
            L.append("    assert MayBeReducible(2, L, info) eq (2 in Set(L));  // family clause rejects p | f unless in L")
    L.append('    printf "SECTION cm: PASS\\n";')
    L.append("end procedure;\n")
    return "\n".join(L)

def _sec_congruence(congs):
    L = ["procedure Test_congruence()"]
    for e in congs:
        eid = e["id"]; strat = e["stratum"]; cong = e["expect"].get("cong")
        L += _build2(e)
        if strat == "fixture-cong-isogenous":
            L.append("    L, info := CongruencePrimes(E1, E2);")
            L.append('    assert info`Source eq "CongruencePrimes";  // %s' % eid)
            L.append('    assert info`Kind eq "AllPrimes";       // isogenous over Q (IsIsogenous)')
            L.append("    L2, info2 := CongruencePrimes(E1, E2 : KnownIsogenous := true);")
            L.append('    assert info2`Kind eq "AllPrimes";')
            L.append('    assert info2`CertificationMethod eq "Supplied";')
        elif strat == "fixture-cong-twist":
            L.append("    L, info := CongruencePrimes(E1, E2);        // %s" % eid)
            L.append('    assert info`Kind eq "Finite";')
            if cong and cong.get("primes") is not None:
                L.append("    assert Set(L) eq %s;" % _mset(cong["primes"]))
            L.append("    L2, info2 := CongruencePrimes(E1, E2 : NormBound := 2, MaxNormBound := 2);")
            L.append('    assert info2`Kind eq "Undecided";     // G = 0 at cap 2; deg-1 short circuits')
            L.append('    assert info2`Kind ne "AllPrimes";     // twist guard: same-j must not certify')
        else:   # fixture-cong-finite
            L.append("    L, info := CongruencePrimes(E1, E2);        // %s" % eid)
            L.append('    assert info`Kind eq "Finite";')
            if cong and cong.get("primes") is not None:
                if eid == "fixture-cong-2" and 2 in cong["primes"]:
                    L.append("    assert 2 in Set(L);  // ell = 2 congruence")
                L.append("    assert Set(L) eq %s;" % _mset(cong["primes"]))
    L.append('    printf "SECTION congruence: PASS\\n";')
    L.append("end procedure;\n")
    return "\n".join(L)

def _sec_fixtures(byid):
    L = ["procedure Test_fixtures()"]
    e = byid.get("fixture-deg1-fldnum")
    if e is not None:
        oE = e["expect"].get("oE")
        L += _build1(e)
        L.append("    assert AbsoluteDegree(BaseRing(E)) eq 1;  // dispatch on absolute degree")
        L.append("    L, info := IsogenyPrimes(E);")
        L.append("    assert info`Exact;                        // branch-1 semantics")
        L.append('    assert info`Kind eq "Finite";')
        if oE is not None:
            L.append("    assert Set(L) eq %s;" % _mset(oE))
    L.append('    printf "SECTION fixtures: PASS\\n";')
    L.append("end procedure;\n")
    return "\n".join(L)

def _dispatch():
    return "\n".join([
        'if section eq "all" or section eq "golden" then Test_golden(); end if;',
        'if section eq "all" or section eq "branch1" then Test_branch1(); end if;',
        'if section eq "all" or section eq "branch2" then Test_branch2(); end if;',
        'if (section eq "all" or section eq "cm") and cmscope ne "0" then Test_cm(); end if;',
        'if section eq "all" or section eq "congruence" then Test_congruence(); end if;',
        'if section eq "all" or section eq "fixtures" then Test_fixtures(); end if;',
        'printf "ALL SELECTED SECTIONS PASS\\n";',
    ])

def emit_test_file(entries, header, path):
    byid = {e["id"]: e for e in entries}
    gates = [e for e in entries if e["id"] in GATE_IDS]
    # Only fixture-cong-* pairs are asserted; diff-congpair is differential-only.
    congs = [e for e in entries if e.get("stratum", "").startswith("fixture-cong")]
    b1, b2, cml = [], [], []
    for e in entries:
        if "curve2" in e:
            continue
        exp = e["expect"]; cm = exp.get("cm") or {}
        if cm.get("is_cm"):
            cml.append(e)
        elif exp.get("deg1"):
            b1.append(e)
        else:
            b2.append(e)
    parts = [PREAMBLE % {"header": header}, "",
             _sec_golden(gates), _sec_branch1(b1), _sec_branch2(b2),
             _sec_cm(cml), _sec_congruence(congs), _sec_fixtures(byid),
             _dispatch(), ""]
    with open(path, "w") as f:
        f.write("\n".join(parts))
    print("wrote %s" % path)

DIFF_HELPERS = r'''
if not assigned out then out := ""; end if;
diffLines := [];

function BoolStr(b)
    return b select "1" else "0";
end function;

function NormPrimes(L)
    S := Sort(SetToSequence(SequenceToSet([ Integers() ! p : p in L ])));
    return Join([ IntegerToString(p) : p in S ], ",");
end function;

procedure EmitIso(~lines, id, L, info)
    if info`IsCM then cmd := IntegerToString(info`CMOrderDiscriminant); else cmd := "0"; end if;
    Append(~lines, id cat ":" cat info`Kind cat ":" cat BoolStr(info`Exact) cat ":" cat NormPrimes(L) cat ":" cat cmd);
end procedure;

procedure EmitCong(~lines, id, L, info)
    Append(~lines, id cat ":" cat info`Kind cat ":" cat BoolStr(info`Exact) cat ":" cat BoolStr(info`Stabilized) cat ":" cat NormPrimes(L) cat ":" cat info`CertificationMethod);
end procedure;
'''

def _diff_entry(e):
    eid = e["id"]
    if "curve2" in e:
        return "\n".join([
            "K := BuildField(%s);" % _mfield(e["field"]),
            "E1 := BuildCurve(K, %s);" % _mainvs(e["curve"]),
            "E2 := BuildCurve(K, %s);" % _mainvs(e["curve2"]),
            "Lc, infoc := CongruencePrimes(E1, E2);",
            'EmitCong(~diffLines, "%s", Lc, infoc);' % eid])
    is_cm = bool((e["expect"].get("cm") or {}).get("is_cm"))
    body = ["K := BuildField(%s);" % _mfield(e["field"]),
            "E := BuildCurve(K, %s);" % _mainvs(e["curve"]),
            "Li, infoi := IsogenyPrimes(E);",
            'EmitIso(~diffLines, "%s", Li, infoi);' % eid]
    if is_cm:
        body = ['if cmscope ne "0" then'] + body + ["end if;"]
    return "\n".join(body)

def emit_differential_driver(entries, header, path):
    parts = [PREAMBLE % {"header": header}, DIFF_HELPERS, ""]
    parts += [_diff_entry(e) for e in entries]
    parts += ["", "Sort(~diffLines);",
              'if out ne "" then',
              '    PrintFile(out, Join(diffLines, "\\n") : Overwrite := true);',
              "else",
              '    for line in diffLines do printf "%o\\n", line; end for;',
              "end if;",
              'printf "DIFFERENTIAL LINES: %o\\n", #diffLines;']
    with open(path, "w") as f:
        f.write("\n".join(parts) + "\n")
    print("wrote %s" % path)

def smoke_subset(entries):
    """Deterministic subset: every fixture-* and sage-doc-* entry, plus the
    first 3 entries of every other stratum. Original order is preserved."""
    keep, cnt = set(), {}
    for e in entries:
        eid = e["id"]
        if eid.startswith("fixture-") or eid.startswith("sage-doc-"):
            keep.add(eid)
        else:
            s = e["stratum"]; cnt[s] = cnt.get(s, 0) + 1
            if cnt[s] <= 3:
                keep.add(eid)
    return [e for e in entries if e["id"] in keep]

def build_header(mode, prov, agg, cap):
    dl = "; ".join("%s:%s" % (i, nm) for (i, nm) in agg) if agg else "none"
    lines = [
        "// mode: %s (per-oracle cap %ds)" % (mode, cap),
        "// sage %s | magma %s | date %s | seed %s" % (
            prov.get("sage_version", "?"), prov.get("magma_version", "?"),
            datetime.date.today().isoformat(), SEED),
        "// LMFDB SQL: %s" % prov.get("sql", "?"),
        "// dropped oracles (entry:oracle): %s" % dl,
    ]
    return "".join(l + "\n" for l in lines)

def _ensure_inputs(d):
    """assemble_inputs() is idempotent (drops and regenerates every non-LMFDB
    entry itself), so this always re-runs it rather than guarding on whether
    fixtures are already present -- a presence guard would skip picking up
    corpus edits (e.g. a fixture pair correction) on a re-run."""
    return assemble_inputs(d)

def main():
    args = sys.argv[1:]
    smoke = "--smoke" in args
    inputs_only = "--inputs-only" in args
    d = json.load(open("tests/corpus_curves.json"))
    d = _ensure_inputs(d)
    if inputs_only:
        json.dump(d, open("tests/corpus_curves.json", "w"), indent=1, default=str)
        print("inputs-only DONE: %d entries" % len(d["entries"]))
        return
    prov = d.get("_provenance", {})
    if smoke:
        entries = smoke_subset(d["entries"])
        print("smoke subset: %d / %d entries (NCPUS=%d)"
              % (len(entries), len(d["entries"]), NCPUS))
        agg = compute_oracles(entries, cap=ORACLE_TIMEOUT)
        header = build_header("smoke", prov, agg, ORACLE_TIMEOUT)
        emit_test_file(entries, header, "/tmp/smoke_test_isogenyprimes.m")
        emit_differential_driver(entries, header, "/tmp/smoke_run_differential.m")
        print("smoke DONE: %d entries, %d dropped" % (len(entries), len(agg)))
    else:
        entries = d["entries"]
        print("full run: %d entries (NCPUS=%d)" % (len(entries), NCPUS))
        agg = compute_oracles(entries, cap=ORACLE_TIMEOUT)
        header = build_header("full", prov, agg, ORACLE_TIMEOUT)
        json.dump(d, open("tests/corpus_curves.json", "w"), indent=1, default=str)
        emit_test_file(entries, header, "tests/test_isogenyprimes.m")
        emit_differential_driver(entries, header, "tests/run_differential.m")
        print("full DONE: %d entries, %d dropped" % (len(entries), len(agg)))

if __name__ == "__main__":
    main()
