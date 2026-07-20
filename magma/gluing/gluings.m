/*
 * Genus2Gluings: the headline intrinsic. Given elliptic curves E1, E2 over Q
 * and a prime n, return every genus-2 curve C over Q whose Jacobian is the
 * (E1 x E2)/graph(psi) gluing along full n-torsion, together with a GluingInfo
 * record of provenance and certificates.
 *
 * Pipeline:
 *   1. require prime n (general n is Tasks 10-11) and base field Q.
 *   2. Congruence prefilter: if GluingModulus is conclusive and n does not
 *      divide it, no gluing exists; return empty with proof "certified".
 *   3. Dispatch:
 *        "Algebraic" -> BHLS closed formulas (n in {2, 3} only);
 *        "Auto"      -> BHLS when n in {2, 3} and the curves are non-isomorphic
 *                       and the BHLS requires hold, else the analytic path;
 *        "Periods"   -> always analytic.
 *   4. Analytic path: GluedPeriodMatrices, RecognizeIgusaClebsch each Jacobian
 *      quotient, reconstruct + pin twist (CurveFromInvariants), CanonicalGluingList.
 *      A near-rational recognition failure doubles the precision and retries
 *      (at most 3 times).
 *
 * Both dispatch paths finish through CanonicalGluingList, so the returned curve
 * list is identical whichever route produced it. proof is "traces-only" until
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

    // Dispatch.
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
    // Phase 1 (recognition, with precision retry): recognize the Igusa-Clebsch
    // invariants of the Jacobian-type quotients. The precision is doubled only
    // while doing so strictly grows the number of recognitions; a plateau with
    // near-rational quotients still unrecognized means those residual quotients
    // are conjugates defined over a number field (Task 8 resolves them) or carry
    // extra automorphisms this reconstruction cannot pin, and further precision
    // will not recover them. The recognition-only loop keeps the (expensive)
    // twist search out of the retry, so n with many quotients stays affordable.
    prec := Precision;
    doublings := 0;
    prevRecognized := -1;
    recIC := []; recPsi := []; products := []; usedPrec := 0;
    while true do
        qs := GluedPeriodMatrices(E1, E2, n : Precision := prec);
        inv0 := qs[1]`invariants;
        usedPrec := fieldPrecision(Parent(inv0[1]));
        recIC := []; recPsi := []; products := []; nearFail := 0;
        for r in qs do
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
            okIC, ICQ := RecognizeIgusaClebsch(r`invariants);
            if okIC then
                Append(~recIC, ICQ); Append(~recPsi, r`psi);
            elif IgusaClebschNearRational(r`invariants) then
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
