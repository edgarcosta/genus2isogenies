/*
 * Shims and small utilities for the gluing package.
 * SiegelReduce wraps the Magma 2.29 builtin whose name is misspelled
 * (To2DUpperHalfSpaceFundamentalDomian); resolved via eval at runtime so the
 * package survives Magma correcting the spelling.
 */

declare verbose Gluing, 2;

intrinsic SiegelReduce(tau::Mtrx) -> Mtrx, Mtrx
{Reduce a 2x2 small period matrix into the Siegel fundamental domain.
Returns the reduced matrix and the integral symplectic transformation.}
    require Nrows(tau) eq 2 and Ncols(tau) eq 2: "tau must be 2x2";
    try
        f := eval "return To2DUpperHalfSpaceFundamentalDomian;";
    catch e
        f := eval "return To2DUpperHalfSpaceFundamentalDomain;";
    end try;
    t, g := f(tau);
    return t, g;
end intrinsic;

intrinsic RecognizeRational(z::FldComElt : Digits := 0) -> BoolElt, FldRatElt
{Attempt to recognize z as a rational number. Heuristic: the imaginary part must
be below 10^(-P/2) and the continued-fraction approximant with denominator at
most 10^(P div 3) must agree with Re(z) to within 10^(-(9*P div 10)), where P is
the working decimal precision (or Digits when given). Returns false, 0 on failure.}
    P := Digits gt 0 select Digits else Floor(Precision(Parent(z)));
    if Abs(Im(z)) gt 10^(-(P div 2)) then return false, Rationals()!0; end if;
    x := Re(z);
    q := BestApproximation(x, 10^(P div 3));
    // Tolerance must sit well below the generic continued-fraction approximation
    // quality (~10^(-2P/3) at this denominator bound, by Dirichlet), else a
    // "generic" irrational (e.g. Pi) gets accepted as often as rejected.
    if Abs(x - q) le 10^(-((9*P) div 10)) * Max(1, Abs(q)) then
        return true, q;
    end if;
    return false, Rationals()!0;
end intrinsic;
