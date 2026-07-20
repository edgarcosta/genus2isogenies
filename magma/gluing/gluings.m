/*
 * Genus2Gluings: the headline intrinsic. Given elliptic curves E1, E2 over Q
 * and a prime n, return every genus-2 curve C over Q whose Jacobian is the
 * (E1 x E2)/graph(psi) gluing along full n-torsion, together with a GluingInfo
 * record of provenance and certificates.
 *
 * Pipeline:
 *   1. require prime n (general n is Tasks 10-11), base field Q, and the
 *      Algorithm parameter's applicability (e.g. "Algebraic" implies n in
 *      {2, 3} and non-isomorphic E1, E2), before anything else runs.
 *   2. Congruence prefilter: if GluingModulus is conclusive and n does not
 *      divide it, no gluing exists; return empty with proof "certified".
 *   3. Dispatch:
 *        "Algebraic" -> BHLS closed formulas (n in {2, 3} only);
 *        "Auto"      -> BHLS when n in {2, 3} and the curves are non-isomorphic
 *                       and the BHLS requires hold, else the analytic path;
 *        "Periods"   -> always analytic.
 *   4. Analytic path: GluedPeriodMatrices (conjugation-filtered), a two-pass
 *      precision sweep (cheap 40-digit look, then RecognizeIgusaClebsch each
 *      surviving Jacobian quotient at working precision), reconstruct + pin
 *      twist (CurveFromInvariants), CanonicalGluingList. A near-rational
 *      recognition failure at working precision doubles it and retries
 *      (at most 3 times).
 *
 * Both dispatch paths finish through CanonicalGluingList; outputs agree only
 * for curves with geometric automorphism group of order 2. Gluings with larger
 * automorphism groups (bielliptic; sextic/quartic twists) are currently produced
 * only by the Algebraic path at n in {2, 3}: CurveFromInvariants detects order > 2
 * on the reconstructed model and drops the quotient with a vprint, since quadratic-
 * twist pinning cannot certify the right twist there. proof is "traces-only" until
 * Task 9 supplies the exact completeness certificate; metadata blocks carry
 * stable_count -1 accordingly.
 */

function gluingInfoFmt()
    return recformat< n : RngIntElt, proof : MonStgElt, blocks : SeqEnum,
                      psis : SeqEnum, products : SeqEnum,
                      precision : RngIntElt, tracebound : RngIntElt >;
end function;

function emptyCurves()
    return [Parent(HyperellipticCurve(PolynomialRing(Rationals()).1^5 + 1)) | ];
end function;

// Working precision of a complex field. A file-local helper so it reaches the
// builtin Precision, which the Genus2Gluings "Precision" parameter shadows
// inside that intrinsic's body.
function fieldPrecision(C)
    return Precision(C);
end function;

// Pass 1 of the two-pass sweep (Genus2Gluings below): a coarse, fixed-
// threshold look at whether a Jacobian-type quotient is plausibly Q-rational,
// cheap enough to run on every conjugation-filtered candidate at the fixed
// low pass-1 precision (40 digits). Mirrors recognize.m's I2-pivot
// normalization (I2 is the only weight dividing 4, 6 and 10) but with a fixed
// gate (1e-10) rather than RecognizeIgusaClebsch/IgusaClebschNearRational's
// precision-scaled one (1e-20 at 40 digits): pass 1's job is to not miss a
// genuine rational quotient to roundoff in a deliberately cheap computation,
// not to make the final call, so it uses the looser threshold.
function LooksRationalIC(IC)
    I2 := IC[1]; I4 := IC[2]; I6 := IC[3]; I10 := IC[4];
    gate := 10^-10;
    if Abs(I10) le gate then return false; end if;
    if Abs(I2^5 / I10) lt gate then return false; end if;
    return Abs(Im(I4 / I2^2)) lt gate and Abs(Im(I6 / I2^3)) lt gate and Abs(Im(I10 / I2^5)) lt gate;
end function;

intrinsic Genus2Gluings(E1::CrvEll, E2::CrvEll, n::RngIntElt
    : Algorithm := "Auto", Precision := false, Proof := "Auto", TraceBound := 1000)
    -> SeqEnum, Rec
{The genus-2 curves over Q gluing E1 and E2 along full n-torsion, for prime n
(general n is lifted in Tasks 10-11), with a GluingInfo record. Algorithm is one
of "Auto" (BHLS closed formulas when n is 2 or 3 and the curves are non-isomorphic,
else analytic periods), "Algebraic" (BHLS, n = 2 or 3 only), or "Periods" (always
analytic). Precision overrides the analytic precision heuristic; TraceBound bounds
the Euler-factor twist certificate. Proof is reserved for the exact certificate of
Task 9. The GluingInfo fields are n, proof ("certified" for the congruence-obstruction
empty answer, else "traces-only"), blocks (one <n, 1, stable_count, analytic_count,
certified> tuple, stable_count -1 until Task 9), psis (a gluing matrix per returned
curve), products (recognized <j1, j2> pairs of the product-type quotients), precision
(digits of the successful analytic pass), and tracebound.}
    require IsPrime(n): "Genus2Gluings currently supports prime n only (general n is Tasks 10-11)";
    require Algorithm in ["Auto", "Algebraic", "Periods"]:
        "Algorithm must be one of \"Auto\", \"Algebraic\", \"Periods\"";
    require Type(BaseRing(E1)) eq FldRat and Type(BaseRing(E2)) eq FldRat:
        "Genus2Gluings currently requires E1 and E2 over Q";

    // Dispatch validation, ahead of the fast path below: an Algorithm choice
    // that does not apply to this (E1, E2, n) must error, even when the pair
    // is certifiably empty (e.g. Algorithm := "Algebraic" at n = 5 must raise
    // the n in {2, 3} require, not silently return an empty list).
    tryAlgebraic := false;
    algebraicStrict := false;   // "Algebraic": let the BHLS requires propagate
    if Algorithm eq "Algebraic" then
        require n in {2, 3}: "the \"Algebraic\" algorithm is only defined for n in {2, 3}";
        require not IsIsomorphic(E1, E2):
            "the \"Algebraic\" algorithm requires non-isomorphic curves";
        tryAlgebraic := true; algebraicStrict := true;
    elif Algorithm eq "Auto" then
        if n in {2, 3} and not IsIsomorphic(E1, E2) then tryAlgebraic := true; end if;
    end if;

    // Congruence-obstruction fast path: a rational n-gluing forces n | a_p(E1) -
    // a_p(E2) at every good p, hence n | GluingModulus. When that modulus is
    // conclusive (not every scanned trace agreed) and n does not divide it, the
    // gluing set is provably empty with no analytic work.
    N, inconclusive := GluingModulus(E1, E2);
    if (not inconclusive) and (N mod n ne 0) then
        info := rec< gluingInfoFmt() | n := n, proof := "certified",
            blocks := [<n, 1, -1, 0, false>], psis := [], products := [],
            precision := 0, tracebound := TraceBound >;
        return emptyCurves(), info;
    end if;

    usedAlgebraic := false;
    cs := emptyCurves();
    if tryAlgebraic then
        if algebraicStrict then
            if n eq 2 then cs := CanonicalGluingList(Genus2Elliptic2(E1, E2));
            else cs := CanonicalGluingList(Genus2Elliptic3(E1, E2)); end if;
            usedAlgebraic := true;
        else
            try
                if n eq 2 then cs := CanonicalGluingList(Genus2Elliptic2(E1, E2));
                else cs := CanonicalGluingList(Genus2Elliptic3(E1, E2)); end if;
                usedAlgebraic := true;
            catch e
                vprintf Gluing: "Genus2Gluings: BHLS unavailable (%o), using periods\n", e`Object;
            end try;
        end if;
    end if;

    if usedAlgebraic then
        info := rec< gluingInfoFmt() | n := n, proof := "traces-only",
            blocks := [<n, 1, -1, #cs, false>], psis := [], products := [],
            precision := 0, tracebound := TraceBound >;
        return cs, info;
    end if;

    // Analytic (periods) path.
    //
    // Phase 1 (recognition), a two-pass precision sweep over the conjugation-
    // filtered candidates GluedPeriodMatrices enumerates:
    //   Pass 1 runs once at a fixed low precision (40 digits, far below the
    //   n-scaled working precision below). Product-type quotients are
    //   recognized (or dropped) right here and never revisited: they never
    //   reach the (expensive) twist search in phase 2 either, so 40 digits is
    //   already enough for their j-invariants. Jacobian-type quotients that
    //   LooksRationalIC flags as plausibly rational are forwarded, by their
    //   psi, to pass 2; one that is visibly complex even at that loose gate is
    //   dropped without ever paying for a full-precision period/theta
    //   computation, which is where the cost of a large candidate pool lives.
    //   Pass 2 recomputes ONLY the forwarded psis, starting at the working
    //   precision (the n/height heuristic, or the caller's Precision
    //   override), and doubles it only while doing so strictly grows the
    //   number of recognitions; a plateau with near-rational quotients still
    //   unrecognized means those residual quotients are conjugates defined
    //   over a number field (Task 8's filter above already removed the
    //   non-equivariant bulk) or carry extra automorphisms this
    //   reconstruction cannot pin, and further precision will not recover
    //   them.
    pass1 := GluedPeriodMatrices(E1, E2, n : Precision := 40);
    products := [];
    survivorPsis := [];
    njacobian := 0;
    for r in pass1 do
        if r`type eq "product" then
            ok1, j1 := RecognizeRational(r`invariants[1]);
            ok2, j2 := RecognizeRational(r`invariants[2]);
            if ok1 and ok2 then
                Append(~products, <j1, j2>);
            else
                vprintf Gluing: "Genus2Gluings: product quotient not over Q, skipping\n";
            end if;
            continue;
        end if;
        njacobian +:= 1;
        if LooksRationalIC(r`invariants) then
            Append(~survivorPsis, r`psi);
        end if;
    end for;
    vprintf Gluing: "Genus2Gluings: pass 1 (precision 40) kept %o/%o jacobian candidate(s) for the full-precision pass\n",
        #survivorPsis, njacobian;

    prec := (Precision cmpeq false) select GluingPrecisionHeuristic(E1, E2, n) else Precision;
    doublings := 0;
    prevRecognized := -1;
    recIC := []; recPsi := []; usedPrec := 40;
    while #survivorPsis gt 0 do
        ws1 := EllipticPeriodBasis(E1, prec);
        ws2 := EllipticPeriodBasis(E2, prec);
        usedPrec := fieldPrecision(Universe(ws1));
        recIC := []; recPsi := []; nearFail := 0;
        for psi in survivorPsis do
            P := GluedBigPeriodMatrix(ws1, ws2, psi, n);
            tau := SmallFromBig(P);
            typ, inv := NumericInvariants(tau);
            if typ eq "product" then
                // Defensive: the product/jacobian gate in NumericInvariants
                // tightens with precision, so a pass-1 jacobian call flipping
                // to product here is not expected (not observed on the
                // corpus); either way it is not a recognizable quotient.
                vprintf Gluing: "Genus2Gluings: pass 2 reclassified a psi as product (unexpected), skipping\n";
                continue;
            end if;
            okIC, ICQ := RecognizeIgusaClebsch(inv);
            if okIC then
                Append(~recIC, ICQ); Append(~recPsi, psi);
            elif IgusaClebschNearRational(inv) then
                nearFail +:= 1;
            end if;
        end for;
        if nearFail eq 0 or doublings ge 3 or #recIC le prevRecognized then
            if nearFail gt 0 then
                vprintf Gluing: "Genus2Gluings: %o near-rational quotient(s) unrecognized (conjugate/unreconstructable), dropping\n", nearFail;
            end if;
            break;
        end if;
        prevRecognized := #recIC;
        prec := 2 * usedPrec;
        doublings +:= 1;
    end while;

    // Phase 2 (reconstruction): recover a Q-model and pin its quadratic twist,
    // once per recognized quotient.
    curves := []; sourcePsis := [];
    for i in [1 .. #recIC] do
        okC, C := CurveFromInvariants(recIC[i], E1, E2 : TraceBound := TraceBound);
        if okC then
            Append(~curves, C); Append(~sourcePsis, recPsi[i]);
        else
            vprintf Gluing: "Genus2Gluings: no twist reproduced L(E1) L(E2), dropping quotient\n";
        end if;
    end for;
    rationalQuotientCount := #curves;

    cs := CanonicalGluingList(curves);
    // A gluing matrix per returned curve. CanonicalGluingList reduces, dedups by
    // isomorphism and sorts, so match each returned curve back to a source psi.
    psisOut := [];
    for c in cs do
        for i in [1 .. #curves] do
            if IsIsomorphic(c, curves[i]) then Append(~psisOut, sourcePsis[i]); break; end if;
        end for;
    end for;
    info := rec< gluingInfoFmt() | n := n, proof := "traces-only",
        blocks := [<n, 1, -1, rationalQuotientCount, false>],
        psis := psisOut, products := products,
        precision := usedPrec, tracebound := TraceBound >;
    return cs, info;
end intrinsic;

intrinsic Genus2Gluings(f1::RngUPolElt, f2::RngUPolElt, n::RngIntElt
    : Algorithm := "Auto", Precision := false, Proof := "Auto", TraceBound := 1000)
    -> SeqEnum, Rec
{Genus2Gluings for the elliptic curves y^2 = f1, y^2 = f2 (f1, f2 univariate
cubics or quartics over Q).}
    return Genus2Gluings(EllipticCurve(f1), EllipticCurve(f2), n
        : Algorithm := Algorithm, Precision := Precision, Proof := Proof, TraceBound := TraceBound);
end intrinsic;

intrinsic Genus2Gluings(s1::MonStgElt, s2::MonStgElt, n::RngIntElt
    : Algorithm := "Auto", Precision := false, Proof := "Auto", TraceBound := 1000)
    -> SeqEnum, Rec
{Genus2Gluings for the elliptic curves with Cremona labels s1, s2.}
    return Genus2Gluings(EllipticCurve(s1), EllipticCurve(s2), n
        : Algorithm := Algorithm, Precision := Precision, Proof := Proof, TraceBound := TraceBound);
end intrinsic;
