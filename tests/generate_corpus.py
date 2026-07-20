"""Corpus generator. assemble_inputs() adds fixture/differential entries to
corpus_curves.json; compute_oracles() (Task 3) adds `expect` blocks and emits
tests/test_isogenyprimes.m and tests/run_differential.m."""
import json
from sage.all import (QQ, ZZ, PolynomialRing, NumberField, EllipticCurve,
                      kronecker_symbol, set_random_seed, prime_range)

SEED = 20260720
MAZUR = [2, 3, 5, 7, 11, 13, 17, 19, 37, 43, 67, 163]

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
    # ell=2 congruence and non-congruent control (Cremona-Freitas style pair
    # is added with its expectation in Task 3 from the 14a/14b trace data).
    es.append(entry("fixture-cong-2", "fixture-cong-finite", QQ,
                    EllipticCurve(QQ, '11a1'), EllipticCurve(QQ, '14a1')))
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

if __name__ == "__main__":
    d = json.load(open("tests/corpus_curves.json"))
    d = assemble_inputs(d)
    json.dump(d, open("tests/corpus_curves.json", "w"), indent=1, default=str)
    print("entries:", len(d["entries"]))
