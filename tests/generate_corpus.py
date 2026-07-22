"""Corpus generator: assembles the test inputs, computes the Sage oracle
expectations, and writes the committed colon-delimited data files that the
hand-written Magma drivers read (tests/test_isogenyprimes.m and
tests/run_differential.m, both via tests/_corpus.m). This module emits DATA
ONLY; the Magma test code is hand-written and lives alongside it in tests/.

Usage:
    sage -python tests/generate_corpus.py            # full run: reads
        tests/corpus_inputs.txt (the fetched LMFDB factors), re-assembles the
        fixture/differential strata, recomputes every Sage oracle, and
        rewrites tests/corpus.txt (minutes: the oE oracle screens each
        (curve, ell) to a certain-negative before any construction, and
        compute_oracles forks across entries). Regression fields are carried
        over from the existing tests/corpus.txt by entry id: they pin engine
        outputs, which a Sage-oracle rerun does not touch.
    sage -python tests/generate_corpus.py --smoke    # deterministic subset,
        writes /tmp/smoke_corpus.txt only, never touches the committed
        corpus (fast; verifies the oracle and writer machinery).
    sage -python tests/generate_corpus.py --add-regressions [out_file]
        # bakes recorded engine outputs into the regression fields of
        tests/corpus.txt. out_file (default .superpowers/sdd/P.out) is in the
        differential .out format written by tests/run_differential.m: a
        curve-pair line is id:kind:exact:stabilized:primes:certmethod, a
        single-curve line id:kind:exact:primes:cmdisc. Text-level rewrite:
        every byte of tests/corpus.txt outside the regression fields and the
        regressions provenance line is preserved.
    sage -python tests/generate_corpus.py --convert-json <corpus_curves.json>
        # one-time conversion from the retired JSON corpus (now in git
        history only): writes tests/corpus_inputs.txt (the fetched LMFDB
        rows) and tests/corpus.txt (every entry with its recorded oracle and
        regression data). A pure format conversion: nothing is recomputed."""
import datetime
import json
import os
import re
import sys
from sage.all import (QQ, ZZ, PolynomialRing, NumberField, EllipticCurve,
                      kronecker_symbol, set_random_seed, prime_range, gcd)
from sage.arith.misc import fundamental_discriminant
from sage.schemes.elliptic_curves.gal_reps_number_field import Billerey_B_l
from sage.parallel.decorate import fork, parallel
from sage.parallel.ncpus import ncpus as _sage_ncpus
try:
    from sage.schemes.elliptic_curves.mod_poly import classical_modular_polynomial
except Exception:      # optional; the modular-polynomial screen degrades to skip
    classical_modular_polynomial = None

SEED = 20260720
MAZUR = [2, 3, 5, 7, 11, 13, 17, 19, 37, 43, 67, 163]
REGRESSION_DEFAULT_OUT = ".superpowers/sdd/P.out"   # --add-regressions default
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
    # strata instead of duplicating them. The lmfdb-qcurve entries (from
    # tests/corpus_inputs.txt) are untouched.
    data["entries"] = [e for e in data["entries"] if e.get("stratum") == "lmfdb-qcurve"]
    es = data["entries"]
    set_random_seed(SEED)
    R = PolynomialRing(QQ, 'x'); x = R.gen()
    # fixture-ex58: Billerey Example 5.8, K = Q(sqrt-3, sqrt-7), conductor 2*O_K.
    # r3, r7 are ARBITRARY roots of x^2+3, x^2+7 in the absolute field (Sage's
    # .roots() order is implementation-defined), so r3*r7 can land on either
    # sign of sqrt(21); only one sign reproduces the paper's curve (the other
    # is an unrelated curve with a huge conductor and nonvanishing B_l, since
    # flipping s21 alone -- unlike flipping s3 or s7 -- is not a field
    # automorphism). Build both signs and keep the one with conductor norm 16,
    # so the fixture is self-verifying instead of trusting root order.
    K58 = NumberField((x**2 + 3, x**2 + 7), names=('s3', 's7')).absolute_field('t')
    r3 = (x**2 + 3).roots(K58, multiplicities=False)[0]
    r7 = (x**2 + 7).roots(K58, multiplicities=False)[0]
    ex58_candidates = []
    for sign in (1, -1):
        r21 = sign * r3 * r7
        a4 = QQ(81)/4 * (69 + 43*r3 + 29*r7 + 17*r21)
        a6 = 162 * (207 - 84*r3 - 54*r7 + 46*r21)
        ex58_candidates.append(EllipticCurve(K58, [0, 0, 0, a4, a6]))
    ex58_matches = [Ec for Ec in ex58_candidates if Ec.conductor().norm() == 16]
    assert len(ex58_matches) == 1, "expected exactly one sqrt(21) sign to give conductor norm 16"
    E58 = ex58_matches[0]
    assert not E58.has_cm()
    # B_l must vanish at every admissible auxiliary prime (r >= 5, r coprime to
    # 6*Disc(K)*Norm(cond E)) -- the fixture's entire point: it forces the
    # engine's B-phase to exhaust its cap and fall through to the R-phase
    # (asserted as `BoundsUsed[2] ne 0` in Test_branch2). Check the first three
    # admissible primes.
    bad58 = ZZ(6) * K58.discriminant() * E58.conductor().norm()
    admissible58 = [r for r in prime_range(1000) if r >= 5 and bad58 % r != 0][:3]
    assert len(admissible58) == 3
    for r in admissible58:
        assert Billerey_B_l(E58, r) == 0, "B_l(%d) != 0 -- wrong ex58 curve" % r
    es.append(entry("fixture-ex58", "fixture-ex58", K58, E58))
    # Re-seed: factoring the huge-conductor rejected candidate and the
    # Billerey_B_l checks above draw from Sage's shared PRNG (current_randstate,
    # e.g. via randomized factoring), which would otherwise shift every
    # ZZ.random_element() draw below (diff-generic-*/diff-x0-*/diff-congpair-*)
    # relative to the pre-fix corpus. Reset to a fresh, deterministic SEED so
    # this fixture's extra computation is invisible to the rest of assembly.
    set_random_seed(SEED)
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
    # fixture-cong-deg1-fldnum: CongruencePrimes on a genuine degree-one FldNum
    # pair (11a1, 11a2, isogenous over Q) must certify AllPrimes via IsIsogenous
    # after transport to Q -- the same absolute-degree dispatch as
    # fixture-deg1-fldnum, exercised on the two-curve entry point.
    es.append(entry("fixture-cong-deg1-fldnum", "fixture-cong-deg1-fldnum", K1,
                    EllipticCurve(K1, [0, -1, 1, -10, -20]),
                    EllipticCurve(K1, [0, -1, 1, -7820, -263580])))
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
    # fixture-h1 (h>1 crash class, found in code review of Task P-7/C-12):
    # EXACT non-minimal Weierstrass models over class-number>1 fields,
    # promoted as regression fixtures for the Reduction-on-non-minimal-model
    # crash fixed in P (e9a9344) and C (7c384a5). Scaling [a1,a2,a3,a4,a6] ->
    # [u a1, u^2 a2, u^3 a3, u^4 a4, u^6 a6] is a genuine Weierstrass change of
    # variable (isomorphic curve, generally non-minimal at primes dividing u),
    # so Sage's own oE/CM oracles are unaffected by it (isogenies_prime_degree/
    # has_cm/cm_discriminant do not require a minimal model): these three flow
    # through the ordinary _compute_expect_for_entry pipeline like any other
    # entry. Only tests/_corpus.m's BuildCurve must reproduce the scaled
    # ainvs
    # exactly, so never call a minimal-model method on these curves. Placed
    # last in assembly (after every seeded/random block above) so nothing here
    # can perturb the shared PRNG state the diff-generic-*/diff-x0-*/
    # diff-congpair-* draws above depend on.
    def _scale_ainvs(ainvs, u):
        a1, a2, a3, a4, a6 = ainvs
        return (u*a1, u**2*a2, u**3*a3, u**4*a4, u**6*a6)
    K5m = NumberField(x**2 + 5, 'w')     # Q(sqrt-5), disc -20, h = 2
    K6m = NumberField(x**2 + 6, 'w')     # Q(sqrt-6), disc -24, h = 2
    assert K5m.class_number() == 2
    assert K6m.class_number() == 2
    # (a) CM j=1728 (order disc -4, i.e. Z[i]) non-minimal at 3.
    Ecm1728 = EllipticCurve(QQ, [0, 0, 0, 1, 0])
    assert Ecm1728.has_cm() and Ecm1728.cm_discriminant() == -4
    a_ainvs = _scale_ainvs(Ecm1728.ainvs(), 3)
    assert a_ainvs == (0, 0, 0, 81, 0)
    es.append(entry("fixture-h1-cm1728", "fixture-h1", K5m,
                    EllipticCurve(K5m, list(a_ainvs))))
    # (b) 11a1 (no CM, conductor 11) non-minimal at the good prime 3.
    E11q = EllipticCurve(QQ, '11a1')
    assert not E11q.has_cm() and E11q.conductor() % 3 != 0
    b_ainvs = _scale_ainvs(E11q.ainvs(), 3)
    assert b_ainvs == (0, -9, 27, -810, -14580)
    es.append(entry("fixture-h1-11a1", "fixture-h1", K5m,
                    EllipticCurve(K5m, list(b_ainvs))))
    # (c) 37a1 (no CM, conductor 37) non-minimal at the good prime 7.
    E37q = EllipticCurve(QQ, '37a1')
    assert not E37q.has_cm() and E37q.conductor() % 7 != 0
    c_ainvs = _scale_ainvs(E37q.ainvs(), 7)
    assert c_ainvs == (0, 0, 343, -2401, 0)
    es.append(entry("fixture-h1-37a1", "fixture-h1", K6m,
                    EllipticCurve(K6m, list(c_ainvs))))
    # fixture-cong-twist-fldnum: degree-two same-j twist BFS gate over
    # Q(sqrt 29). E = 11a1 base-changed (E11 above); E' = its quadratic twist by
    # -1 (all-rational model [0,4,0,-160,1264], same j, non-isomorphic). The
    # first prime good for both curves whose Frobenius charpolies differ has
    # norm 7, so CongruencePrimes at cap 6 (below 7, but above the two good
    # norm-5 primes) reaches G = 0 and the degree-two certification BFS runs; it
    # must REFUSE to certify the same-j-but-non-isomorphic node (Undecided,
    # never AllPrimes), the case a degree-one short circuit ([_,0,0]) and a
    # j-only matcher (AllPrimes) both get wrong. Placed last (like fixture-h1)
    # so its conductor/factoring calls cannot perturb the shared PRNG the
    # diff-* draws above depend on.
    E11tw = E11.quadratic_twist(K29(-1))
    assert list(E11tw.ainvs()) == [0, 4, 0, -160, 1264]
    assert not K29(-1).is_square()
    assert E11.j_invariant() == E11tw.j_invariant()
    assert not E11.is_isomorphic(E11tw)
    _Tw = PolynomialRing(ZZ, 'Tw').gen()
    def _twist_fcp(E, P):
        return _Tw**2 - ZZ(E.reduction(P).trace_of_frobenius()) * _Tw + ZZ(P.norm())
    _tw_places = [P for p in prime_range(2, 50) for P in K29.primes_above(p)]
    _tw_common_good = [P for P in _tw_places
                       if E11.conductor().valuation(P) == 0
                       and E11tw.conductor().valuation(P) == 0]
    _tw_sep = [P for P in _tw_common_good
               if _twist_fcp(E11, P) != _twist_fcp(E11tw, P)]
    assert min(ZZ(P.norm()) for P in _tw_sep) == 7   # first separating norm
    assert all(_twist_fcp(E11, P) == _twist_fcp(E11tw, P)
               for P in _tw_common_good if P.norm() <= 6)   # cap 6 forces G = 0
    es.append(entry("fixture-cong-twist-fldnum", "fixture-cong-twist-fldnum",
                    K29, E11, E11tw))
    return data

# ===========================================================================
# Oracle computation.
#
# The oracles below reconstruct each entry's K and E in Sage the same way the
# Magma BuildField/BuildCurve constructors in tests/_corpus.m do. Sage and
# Magma may pick
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
    silently. A Sage exception in the construction arm records the same
    (curve, ell) drop, never a silent negative. `stats` is filled with the
    per-screen hit counts for logging."""
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
                # A Sage arithmetic failure is NO DATA for this ell, not
                # negative evidence (spec: drops are never silent). Recorded by
                # the caller exactly like a timeout drop, never a silent "out".
                return "error_drop"
        r = capped("oE:ell=%d" % elli, per_ell, dropped)
        if r is None:                 # timeout: capped() already recorded the drop
            n_drop += 1
        elif r == "error_drop":       # exception inside the oracle: same drop bookkeeping
            dropped.append("oE:ell=%d" % elli)
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
    it). Returns the aggregated [(id, oracle)] drop list; the corpus writer
    re-derives the same aggregate for the header provenance.
    Each entry's expect block is a pure function of the entry, so results are
    reattached BY ID in the original entry order: the written corpus is
    byte-identical regardless of the order the forks finish. A
    whole-entry timeout (or crash) is recorded as an "entry_timeout" drop rather
    than crashing compute_oracles or silently omitting the entry: the entry
    keeps its id and gets a minimal expect block, so the absent oracles are
    written as empty fields and the drivers skip its asserts cleanly instead
    of asserting something wrong."""
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
            # drivers' branch1-vs-branch2 classification still matches the
            # entry's real shape; every other oracle is left absent, which
            # the corpus writer records as empty fields and the Magma
            # drivers treat as "oracle dropped".
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

# ===========================================================================
# The committed data files. tests/corpus.txt is the single data file both
# Magma drivers read (schema documented in its header comment, parsed by
# tests/_corpus.m); tests/corpus_inputs.txt is the fetched LMFDB slice a full
# run re-expands into the lmfdb-qcurve stratum. This module is the only
# writer of either file. House format (data/README.md): colon-delimited
# fields, no spaces, bracketed lists; a colon never appears at bracket depth
# 0 inside a field, so consumers split lines on top-level colons.
# ===========================================================================

CORPUS_PATH = "tests/corpus.txt"
INPUTS_PATH = "tests/corpus_inputs.txt"

def split_top(s, sep=":"):
    """Split s on sep at bracket depth 0, preserving empty fields (Python's
    str.split would break bracketed lists; Magma's Split also drops empty
    fields, so tests/_corpus.m carries this same loop)."""
    parts, cur, depth = [], [], 0
    for ch in s:
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
        if ch == sep and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return parts

def _bool01(b):
    return "1" if b else "0"

def _blist(items):
    return "[" + ",".join(items) + "]"

def _vec(v):
    return _blist([str(c) for c in v])

def _curve_rows(model):
    if model["model"] == "AB":
        return [[0], [0], [0], model["A"], model["B"]]
    return model["ainvs"]

def _curve_field(model):
    return _blist([_vec(r) for r in _curve_rows(model)])

def _ints(v):
    """Bracketed integer list; None (dropped or unavailable oracle) is the
    empty field, distinct from [] (computed and empty)."""
    if v is None:
        return ""
    return _blist([str(int(x)) for x in v])

def _cm_field(cm):
    if cm is None:
        return ""
    if not cm.get("is_cm"):
        return "[0]"
    return _blist(["1", str(cm["order_disc"]), str(cm["fund_disc"]),
                   str(cm["f"]), _bool01(cm["in_base_field"])])

def _golden_field(g):
    if g is None:
        return ""
    rows = []
    for tag in ("split", "inert"):
        r = g[tag]
        rows.append(_blist([str(r["p"]), str(r["norm"])]
                           + [str(c) for c in r["charpoly"]]))
    return _blist(rows)

def _cong_fields(c):
    if c is None:
        return ["", "", "", ""]
    return [c["kind"], _bool01(c["stabilized"]), str(c["bound"]),
            _ints(c["primes"])]

def _reg_single(r):
    if r is None:
        return ["", "", "", ""]
    return [r["kind"], _bool01(r["exact"]), _ints(r["primes"]),
            str(r["cmdisc"])]

def _reg_pair(r):
    if r is None:
        return ["", "", "", "", ""]
    return [r["kind"], _bool01(r["exact"]), _bool01(r["stabilized"]),
            _ints(r["primes"]), r["certmethod"]]

def entry_line(e, reg_override=None):
    """One corpus.txt data line. reg_override, when given, is the raw list of
    regression field strings harvested from a previous corpus.txt (full runs
    carry pins forward by id; a fresh entry has none)."""
    exp = e.get("expect") or {}
    deg1 = exp.get("deg1")
    if deg1 is None:
        # entry_timeout fallback; matches the parser's field-degree cross-check
        deg1 = len(e["field"]) == 2
    f = [e["id"], e["stratum"], _vec(e["field"]), _curve_field(e["curve"])]
    if "curve2" in e:
        f.append(_curve_field(e["curve2"]))
        f.append(_bool01(deg1))
        f += _cong_fields(exp.get("cong"))
        f += _cong_fields(exp.get("cong_lowbound"))
        f.append(_blist(list(exp.get("dropped") or [])))
        f += reg_override if reg_override is not None else _reg_pair(e.get("regression"))
        assert len(f) == 20, (e["id"], len(f))
    else:
        f.append(_bool01(deg1))
        f.append(_cm_field(exp.get("cm")))
        f.append(_ints(exp.get("oE")))
        f.append(_ints(exp.get("soft_reducible")))
        f.append(_golden_field(exp.get("golden")))
        f.append(_blist(list(exp.get("dropped") or [])))
        f += reg_override if reg_override is not None else _reg_single(e.get("regression"))
        assert len(f) == 14, (e["id"], len(f))
    line = ":".join(f)
    assert split_top(line) == f, "field with a top-level colon: %r" % line
    return line

def _agg_drops(entries):
    agg = []
    for e in entries:
        expect = e.get("expect") or {}
        agg.extend((e["id"], nm) for nm in (expect.get("dropped") or []))
    return agg

_SCHEMA_LINES = """\
# Test corpus for the isogeny-primes engine (magma/IsogenyPrimes.m). Read by
# tests/_corpus.m for tests/test_isogenyprimes.m and tests/run_differential.m;
# written by tests/generate_corpus.py (do not edit data lines by hand).
#
# One entry per line. Fields are colon-delimited with no spaces; lists are
# bracketed as in data/README.md. A colon never appears at bracket depth 0
# inside a field, so lines split on top-level colons (dropped-oracle labels
# keep their colons inside the bracketed list). An empty field means the
# value is absent (oracle dropped, not applicable, or not pinned), which is
# distinct from [] (computed and empty).
#
# Single-curve entries, 14 fields:
#   id:stratum:field:ainvs:deg1:cm:oE:soft:golden:dropped:regkind:regexact:regprimes:regcmdisc
#   field    ascending coefficients of the defining polynomial of K over Q;
#            [0,1] is Q itself, any other degree-1 polynomial builds the
#            genuine degree-one FldNum (RationalsAsNumberField).
#   ainvs    [a1,a2,a3,a4,a6], each a power-basis coordinate vector over Q
#            (LMFDB short-Weierstrass factors appear as [[0],[0],[0],A,B]).
#   deg1     1 iff K has absolute degree 1; redundant with #field, kept as a
#            hand-edit tripwire (the parser cross-checks them).
#   cm       Sage geometric-CM oracle: [0] no CM; [1,DO,DF,f,ib] CM with
#            order discriminant DO = f^2*DF, fundamental discriminant DF,
#            conductor f, and ib = 1 iff sqrt(DF) lies in K.
#   oE       Sage construction oracle O(E): the sorted ells with
#            isogenies_prime_degree(ell) nonempty (the Mazur list over
#            degree-one fields, ell < 101 otherwise).
#   soft     Sage reducible_primes(), print-only tightness comparison
#            (absent on degree-one and CM entries).
#   golden   Sage Frobenius oracle for the two gate entries, one
#            [p,Norm(q),c0,c1,c2] row per prime (split then inert): the
#            charpoly at a prime q above p has ascending coefficients
#            c0,c1,c2.
#   dropped  per-entry dropped-oracle labels (fork timeout or error), e.g.
#            oE:ell=97; any oE:* drop makes the oE oracle a lower bound, so
#            equality asserts downgrade to containment.
#   reg*     recorded engine output of IsogenyPrimes (regression pin) in the
#            differential line vocabulary: Kind, Exact (0/1), sorted primes,
#            CM order discriminant (0 when not CM). Empty regkind = this
#            entry's engine output is not pinned.
#
# Curve-pair entries, 20 fields:
#   id:stratum:field:ainvs:ainvs2:deg1:congkind:congstab:congbound:congprimes:lowkind:lowstab:lowbound:lowprimes:dropped:regkind:regexact:regstab:regprimes:regcert
#   cong*    Sage trace-gcd congruence oracle at the default escalation
#            (norm bound 1000 doubling to 8000): kind Finite or ZeroAtCap,
#            stabilized (0/1), final bound, prime divisors of the final gcd
#            (empty congprimes with kind ZeroAtCap: the gcd was 0 at cap).
#            Empty congkind = no congruence oracle (differential-only pair).
#   low*     the same oracle rerun capped at norm bound 2 (fixture-cong-twist
#            only): pins the zero-sentinel expectation for the low-cap
#            engine call.
#   reg*     recorded engine output of CongruencePrimes: Kind, Exact,
#            Stabilized, sorted primes, CertificationMethod.
#""".splitlines()

def _regressions_line(count, magma_ver, baked):
    if not count:
        return ("# regressions: none recorded "
                "(run tests/run_differential.m, then --add-regressions)")
    return ("# regressions: %d recorded engine outputs "
            "(differential line format, magma %s), baked %s"
            % (count, magma_ver or "unknown", baked or "unknown"))

def write_corpus(entries, prov, path=CORPUS_PATH, reg_overrides=None):
    reg_overrides = reg_overrides or {}
    agg = _agg_drops(entries)
    dl = "; ".join("%s:%s" % (i, nm) for (i, nm) in agg) if agg else "none"
    npins = 0
    body = []
    for e in entries:
        line = entry_line(e, reg_override=reg_overrides.get(e["id"]))
        fields = split_top(line)
        if (fields[15] if len(fields) == 20 else fields[10]) != "":
            npins += 1
        body.append(line)
    head = list(_SCHEMA_LINES)
    head.append("# provenance: oracles sage %s, seed %s, per-oracle cap %ds, computed %s"
                % (prov.get("sage_version", "?"), SEED, ORACLE_TIMEOUT,
                   prov.get("oracle_date", "?")))
    head.append(_regressions_line(npins, prov.get("reg_magma"),
                                  prov.get("reg_baked")))
    head.append("# lmfdb: fetched %s, source %s"
                % (prov.get("fetched", "?"), prov.get("source", "?")))
    head.append("# lmfdb sql: %s" % prov.get("sql", "?"))
    head.append("# dropped oracles (entry:oracle): %s" % dl)
    with open(path, "w") as fh:
        fh.write("\n".join(head + body) + "\n")
    print("wrote %s: %d entries, %d pinned" % (path, len(body), npins))

# --- tests/corpus_inputs.txt -------------------------------------------------

_INPUTS_SCHEMA = """\
# LMFDB inputs for the isogeny-primes test corpus: the elliptic factors, over
# their quadratic splitting fields, of the LMFDB genus 2 curves whose
# Jacobian is geometrically Q x Q but splits only over a quadratic field.
# Read and re-expanded into the lmfdb-qcurve stratum by
# tests/generate_corpus.py; regenerated by the corpus fetcher.
# One factor per line, colon-delimited, no spaces:
#   label:factor:spl_fod_label:factor_label:field:curve
#   label          g2c curve label
#   factor         index of this factor in the row's spl_facs_coeffs
#   spl_fod_label  LMFDB label of the quadratic splitting field
#   factor_label   LMFDB label of this elliptic factor over that field
#                  (empty when the LMFDB has none)
#   field          spl_fod_coeffs, ascending
#   curve          this factor's spl_facs_coeffs pair [A,B]: the curve
#                  y^2 = x^3 + A*x + B with A, B power-basis coordinate
#                  vectors over Q""".splitlines()

def _tokens(s):
    """The comma-separated top-level tokens of one bracketed list ('[]' gives
    the empty list). Tokens stay textual: they round-trip to the writer and
    coerce through QQ() on the Sage side."""
    assert s.startswith("[") and s.endswith("]"), s
    inner = s[1:-1]
    return split_top(inner, ",") if inner else []

def write_inputs(rows, fetched, source, sql, path=INPUTS_PATH):
    head = list(_INPUTS_SCHEMA)
    head.append("# fetched: %s, source %s" % (fetched, source))
    head.append("# sql: %s" % sql)
    body = []
    for r in rows:
        line = ":".join(r)
        assert split_top(line) == list(r), "field with a top-level colon: %r" % line
        body.append(line)
    with open(path, "w") as fh:
        fh.write("\n".join(head + body) + "\n")
    print("wrote %s: %d factor rows" % (path, len(body)))

def read_inputs(path=INPUTS_PATH):
    prov = {}
    entries = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                m = re.match(r"# fetched: (\S+), source (.*)$", line)
                if m:
                    prov["fetched"], prov["source"] = m.group(1), m.group(2)
                m = re.match(r"# sql: (.*)$", line)
                if m:
                    prov["sql"] = m.group(1)
                continue
            f = split_top(line)
            assert len(f) == 6, "inputs line with %d fields: %r" % (len(f), line)
            label, idx, fod, faclabel, field, curve = f
            AB = _tokens(curve)
            assert len(AB) == 2, "curve field is not an [A,B] pair: %r" % line
            entries.append({
                "id": "lmfdb-%s-%s" % (label, idx),
                "stratum": "lmfdb-qcurve",
                "lmfdb_factor_label": faclabel or None,
                "spl_fod_label": fod,
                "field": _tokens(field),
                "curve": {"model": "AB", "A": _tokens(AB[0]), "B": _tokens(AB[1])},
            })
    assert entries, "no data lines in %s" % path
    return prov, entries

# --- regression baking (text-level rewrite of tests/corpus.txt) --------------

def _reg_span(nfields):
    """(start, stop) slice of the regression fields for a data line with
    nfields fields (14 = single curve, 20 = pair)."""
    assert nfields in (14, 20), nfields
    return (15, 20) if nfields == 20 else (10, 14)

def harvest_regressions(path=CORPUS_PATH):
    """Regression fields by id from an existing corpus.txt, plus the parsed
    provenance of its regressions header line. Pins record engine outputs, so
    a Sage-oracle rerun carries them forward untouched."""
    if not os.path.exists(path):
        return {}, (None, None)
    over = {}
    reg_magma = reg_baked = None
    for line in open(path).read().splitlines():
        if line.startswith("# regressions:"):
            m = re.search(r"magma ([^)]+)\), baked (\S+)", line)
            if m:
                reg_magma, reg_baked = m.group(1), m.group(2)
            continue
        if not line or line.startswith("#"):
            continue
        f = split_top(line)
        a, b = _reg_span(len(f))
        if f[a] != "":
            over[f[0]] = f[a:b]
    return over, (reg_magma, reg_baked)

def add_regressions(out_path):
    lines = open(CORPUS_PATH).read().splitlines()
    ispair, lineno = {}, {}
    for i, line in enumerate(lines):
        if not line or line.startswith("#"):
            continue
        f = split_top(line)
        _reg_span(len(f))          # validates the field count
        ispair[f[0]] = (len(f) == 20)
        lineno[f[0]] = i
    n = 0
    with open(out_path) as fh:
        for raw in fh:
            raw = raw.rstrip("\n")
            if not raw:
                continue
            parts = raw.split(":")
            eid = parts[0]
            assert eid in ispair, (
                "add-regressions: id %r not found in corpus (%r)" % (eid, raw))
            f = split_top(lines[lineno[eid]])
            if ispair[eid]:
                assert len(parts) == 6, (
                    "add-regressions: %r is a pair entry but line has %d "
                    "fields, expected 6: %r" % (eid, len(parts), raw))
                _, kind, exact, stabilized, primes, certmethod = parts
                f[15:20] = [kind, _bool01(int(exact)), _bool01(int(stabilized)),
                            _ints([int(p) for p in primes.split(",")] if primes else []),
                            certmethod]
            else:
                assert len(parts) == 5, (
                    "add-regressions: %r is a single-curve entry but line has "
                    "%d fields, expected 5: %r" % (eid, len(parts), raw))
                _, kind, exact, primes, cmdisc = parts
                f[10:14] = [kind, _bool01(int(exact)),
                            _ints([int(p) for p in primes.split(",")] if primes else []),
                            str(int(cmdisc))]
            lines[lineno[eid]] = ":".join(f)
            n += 1
    npins = 0
    for line in lines:
        if not line or line.startswith("#"):
            continue
        f = split_top(line)
        a, _b = _reg_span(len(f))
        if f[a] != "":
            npins += 1
    ver = None
    for i, line in enumerate(lines):
        if line.startswith("# regressions:"):
            m = re.search(r"magma ([^)]+)\)", line)
            if m:
                ver = m.group(1)
            lines[i] = _regressions_line(npins, ver,
                                         datetime.date.today().isoformat())
            break
    with open(CORPUS_PATH, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("add-regressions DONE: %d records applied from %s; %d entries pinned"
          % (n, out_path, npins))

# --- one-time conversion from the retired JSON corpus ------------------------

def _fund_disc(D):
    """Fundamental discriminant of Q(sqrt D) for a nonsquare integer D: the
    squarefree part of D, times 4 unless that is 1 mod 4."""
    s = -1 if D < 0 else 1
    n = abs(D)
    d = 2
    while d * d <= n:
        while n % (d * d) == 0:
            n //= d * d
        d += 1
    d0 = s * n
    return d0 if d0 % 4 == 1 else 4 * d0

def _quad_field_label(coeffs):
    """LMFDB label (2.r1.|disc|.1) of the quadratic field with the given
    monic defining polynomial."""
    c = [int(x) for x in coeffs]
    assert len(c) == 3 and c[2] == 1, coeffs
    d0 = _fund_disc(c[1] * c[1] - 4 * c[0])
    return "2.%d.%d.1" % (2 if d0 > 0 else 0, abs(d0))

def convert_json(json_path):
    """One-time conversion of the retired JSON corpus to the committed data
    files. Pure format conversion: every oracle expectation and regression
    pin is carried over verbatim, nothing is recomputed. The JSON did not
    store spl_fod_label, so it is derived from the field discriminant and
    cross-checked against the factor-label prefix wherever the LMFDB
    supplied one."""
    d = json.load(open(json_path))
    entries = d["entries"]
    jprov = d.get("_provenance", {})
    rows = []
    for e in entries:
        if e.get("stratum") != "lmfdb-qcurve":
            continue
        label, idx = e["id"][len("lmfdb-"):].rsplit("-", 1)
        fod = _quad_field_label(e["field"])
        faclabel = e.get("lmfdb_factor_label") or ""
        if faclabel:
            assert faclabel.split("-")[0] == fod, (e["id"], faclabel, fod)
        rows.append((label, idx, fod, faclabel, _vec(e["field"]),
                     _blist([_vec(e["curve"]["A"]), _vec(e["curve"]["B"])])))
    write_inputs(rows, jprov.get("date", "?"), jprov.get("source", "?"),
                 jprov.get("sql", "?"))
    prov = {
        "sage_version": jprov.get("sage_version", "?"),
        # The committed oracle data and regression pins were last recomputed
        # and re-baked on 2026-07-22 (the emission date of the generated
        # files this conversion supersedes).
        "oracle_date": "2026-07-22",
        "reg_magma": jprov.get("magma_version", "unknown"),
        "reg_baked": "2026-07-22",
        "fetched": jprov.get("date", "?"),
        "source": jprov.get("source", "?"),
        "sql": jprov.get("sql", "?"),
    }
    write_corpus(entries, prov)

# --- entry points ------------------------------------------------------------

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

def _sage_version():
    try:
        from sage.version import version
        return version
    except Exception:
        return "?"

def main():
    args = sys.argv[1:]
    if "--convert-json" in args:
        idx = args.index("--convert-json")
        assert idx + 1 < len(args), "--convert-json needs the JSON path"
        convert_json(args[idx + 1])
        return
    if "--add-regressions" in args:
        idx = args.index("--add-regressions")
        out_path = REGRESSION_DEFAULT_OUT
        if idx + 1 < len(args) and not args[idx + 1].startswith("--"):
            out_path = args[idx + 1]
        add_regressions(out_path)
        return
    smoke = "--smoke" in args
    iprov, lm = read_inputs()
    d = assemble_inputs({"entries": lm})
    prov = {
        "sage_version": _sage_version(),
        "oracle_date": datetime.date.today().isoformat(),
        "fetched": iprov.get("fetched", "?"),
        "source": iprov.get("source", "?"),
        "sql": iprov.get("sql", "?"),
    }
    if smoke:
        entries = smoke_subset(d["entries"])
        print("smoke subset: %d / %d entries (NCPUS=%d)"
              % (len(entries), len(d["entries"]), NCPUS))
        compute_oracles(entries, cap=ORACLE_TIMEOUT)
        write_corpus(entries, prov, path="/tmp/smoke_corpus.txt")
        print("smoke DONE")
    else:
        entries = d["entries"]
        print("full run: %d entries (NCPUS=%d)" % (len(entries), NCPUS))
        compute_oracles(entries, cap=ORACLE_TIMEOUT)
        over, (reg_magma, reg_baked) = harvest_regressions()
        prov["reg_magma"], prov["reg_baked"] = reg_magma, reg_baked
        write_corpus(entries, prov, reg_overrides=over)
        print("full DONE: %d entries" % len(entries))

if __name__ == "__main__":
    main()
