/*
 * Congruence-modulus prefilter: a necessary condition on n for a rational
 * gluing of E1 and E2 along full n-torsion to exist, cheap to test before
 * running the analytic or algebraic gluing machinery.
 *
 * Also home to GeometricallyIsogenous, the exact scope test for the paper's
 * Algorithm 3.4 (which takes geometrically NONisogenous pairs): the certificate
 * layer (gluings.m) abstains whenever it returns true.
 */

// Exact CM detection over Q: E has complex multiplication iff j(E) is one of the
// 13 rational CM j-invariants, i.e. a root of HilbertClassPolynomial(D) for one of
// the 13 imaginary quadratic orders D of class number one. A finite table lookup,
// not an analytic decision; returns the discriminant on a hit.
function rationalCMDiscriminant(j)
    for D in [-3, -4, -7, -8, -11, -12, -16, -19, -27, -28, -43, -67, -163] do
        if Evaluate(HilbertClassPolynomial(D), j) eq 0 then return true, D; end if;
    end for;
    return false, 0;
end function;

intrinsic GeometricallyIsogenous(E1::CrvEll, E2::CrvEll) -> BoolElt
{True iff E1 and E2 are isogenous over Qbar, i.e. Hom_Qbar(E1, E2) is nonzero. The
test is exact, not heuristic. Equal j-invariants: geometrically isomorphic, true.
Otherwise CM is decided by comparing j against the 13 rational CM j-invariants
(HilbertClassPolynomial of the class-number-one discriminants): two CM curves are
geometrically isogenous iff their CM fields agree (product of discriminants a
square), and a CM curve is never geometrically isogenous to a non-CM curve (an
isogeny preserves End tensor Q). Two non-CM curves are geometrically isogenous iff
some member of the Q-isogeny class of E1 has the j-invariant of E2: the geometric
Hom generator of a non-CM pair is rational after a quadratic twist, so such a
member exists exactly when Hom_Qbar(E1, E2) is nonzero (the IsQuadraticTwist call
is a defensive consistency check only). E1 and E2 must be defined over Q.}
    require Type(BaseRing(E1)) eq FldRat and Type(BaseRing(E2)) eq FldRat:
        "GeometricallyIsogenous requires E1 and E2 defined over Q";
    j1 := jInvariant(E1); j2 := jInvariant(E2);
    if j1 eq j2 then return true; end if;
    cm1, D1 := rationalCMDiscriminant(j1);
    cm2, D2 := rationalCMDiscriminant(j2);
    if cm1 or cm2 then
        return cm1 and cm2 and IsSquare(D1 * D2);
    end if;
    for F in IsogenousCurves(E1) do
        if jInvariant(F) eq j2 then
            ok, _ := IsQuadraticTwist(F, E2);
            if not ok then
                vprintf Gluing: "GeometricallyIsogenous: same-j member of the isogeny class is not a quadratic twist of E2 (unexpected for non-CM); the j-match already decides\n";
            end if;
            return true;
        end if;
    end for;
    return false;
end intrinsic;

intrinsic GluingModulus(E1::CrvEll, E2::CrvEll : Bound := 500) -> RngIntElt, BoolElt
{N = gcd over primes p <= Bound of good reduction for both curves of
a_p(E1) - a_p(E2). Any rational gluing along full n-torsion requires n | N
(Algorithm split-non-square, Steps 1-2). Second return value is true when every
scanned prime gave equal traces: the curves are then plausibly isogenous and the
returned N = 0 carries no information.}
    K := BaseRing(E1);
    require BaseRing(E2) cmpeq K: "curves must share a base field";
    N := 0;
    if Type(K) eq FldRat then
        bad := Set(BadPrimes(E1)) join Set(BadPrimes(E2));
        for p in PrimesUpTo(Bound) do
            if p in bad then continue; end if;
            N := GCD(N, TraceOfFrobenius(E1, p) - TraceOfFrobenius(E2, p));
            if N eq 1 then return 1, false; end if;
        end for;
    else
        c1 := Conductor(E1); c2 := Conductor(E2);
        OK := MaximalOrder(K);
        for p in PrimesUpTo(Bound) do
            for P in [f[1] : f in Factorization(p * OK)] do
                if Valuation(c1, P) ne 0 or Valuation(c2, P) ne 0 then continue; end if;
                N := GCD(N, TraceOfFrobenius(E1, P) - TraceOfFrobenius(E2, P));
                if N eq 1 then return 1, false; end if;
            end for;
        end for;
    end if;
    return Abs(N), N eq 0;
end intrinsic;
