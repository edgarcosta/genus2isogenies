/*
 * Congruence-modulus prefilter: a necessary condition on n for a rational
 * gluing of E1 and E2 along full n-torsion to exist, cheap to test before
 * running the analytic or algebraic gluing machinery.
 */

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
