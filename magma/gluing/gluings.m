/*
 * Genus2Gluings: the headline intrinsic. Given elliptic curves E1, E2 over Q
 * and any n >= 2, return every genus-2 curve C over Q whose Jacobian is the
 * (E1 x E2)/graph(psi) gluing along full n-torsion, together with a GluingInfo
 * record of provenance and certificates. A composite n is composed from its
 * prime-power blocks by CRT (gluingCompositeCRT); a prime power n = ell^e runs
 * the pipeline below directly.
 *
 * Pipeline:
 *   1. require n >= 2, base field Q, and the Algorithm parameter's applicability
 *      (e.g. "Algebraic" implies n in {2, 3} and non-isomorphic E1, E2), before
 *      anything else runs. Composite n dispatches to gluingCompositeCRT here.
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
 * twist pinning cannot certify the right twist there.
 *
 * The exact completeness certificate (exact.m) closes the loop: GaloisStableGluings
 * counts the Galois-stable anti-symplectic gluings exactly (as M/-M-orbit quotients),
 * and GluingCertificateBlock compares that against the analytic count of distinct
 * recognized rational quotients (#Seqset(recIC) + #Seqset(products) from
 * gluingAnalyticEnumerate, which counts every recognized quotient including the
 * bielliptic/Mestre ones this path cannot emit and the product-type ones, one per
 * moduli class). Under Proof := "Auto" (prime n <= 13) or true, a match certifies
 * the block and a disagreement is a hard error; the proof field is "certified" iff
 * every block certified.
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

// The two-pass precision sweep, shared by the analytic dispatch (which
// reconstructs curves from the result) and the completeness certificate (which
// only needs the count). Enumerates the conjugation-filtered anti-symplectic
// quotients GluedPeriodMatrices returns and recognizes the rational ones: pass 1
// at a fixed 40 digits classifies products (recognized here) and forwards
// plausibly-rational jacobian psi; pass 2 recomputes those at working precision,
// doubling only while that strictly grows the recognition count. Returns the
// recognized rational jacobian Igusa-Clebsch list recIC and their source psi
// recPsi (per psi: both psi and -psi appear), the recognized rational product
// j-pairs products (per psi), and the working precision used.
//
// The certificate's analytic count is #Seqset(recIC) + #Seqset(products): DISTINCT
// recognized quotients (deduplicated by moduli), the same unit as the exact side's
// M/-M-orbit count. It counts a recognized jacobian quotient whether or not its
// curve model survives twist pinning in phase 2 (a recognized rational quotient is
// Galois-stable regardless of the bielliptic/Mestre reconstruction obstructions),
// but deduplicates so that two psi landing on the same rational moduli (psi and
// -psi always; also look-alike graphs when the target has extra automorphisms)
// count once, matching how the exact side collapses graphs to quotients.
function gluingAnalyticEnumerate(E1, E2, n, prec)
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
        if LooksRationalIC(r`invariants) then Append(~survivorPsis, r`psi); end if;
    end for;
    vprintf Gluing: "Genus2Gluings: pass 1 (precision 40) kept %o/%o jacobian candidate(s) for the full-precision pass\n",
        #survivorPsis, njacobian;

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
    return recIC, recPsi, products, usedPrec;
end function;

// Does E admit a Q-rational ell-isogeny (a Galois-stable order-ell subgroup)?
// A product ell-gluing is E1' x E2' with E1' ell-isogenous to E1 and E2' to E2
// over Q, so if either curve has no rational ell-isogeny there is NO Galois-stable
// product gluing of the pair. A recognized rational product is then a look-alike
// (a non-Galois-stable graph whose product quotient happens to share moduli with
// a Q-rational product that does not descend through this graph, as at j = 0 with
// an irreducible mod-ell image) and must be dropped from the certificate count so
// it matches the exact Galois-stable count. When both curves do admit a rational
// ell-isogeny the recognized products are kept; any residual over-count then
// surfaces as a loud certificate mismatch rather than a silent false certification.
function admitsRationalIsogeny(E, ell)
    if ell eq 2 then
        return #Roots(DivisionPolynomial(E, 2)) gt 0;   // rational 2-torsion point
    end if;
    // odd ell: a rational point on X_0(ell) over j(E), i.e. a rational root of
    // Phi_ell(j(E), Y) (the ell-th classical modular polynomial). Phi_ell is
    // integral, so change to Q before substituting the (rational) j-invariant.
    Phi := ChangeRing(ClassicalModularPolynomial(ell), Rationals());
    return #Roots(UnivariatePolynomial(Evaluate(Phi, 1, jInvariant(E)))) gt 0;
end function;

// admitsRationalIsogeny, but for a level past ClassicalModularPolynomial's database
// (ell > 59): that lookup errors rather than returning false, so treat "unavailable"
// as "cannot rule out" (gate passes, product kept) instead of crashing.
function admitsRationalIsogenyOrUnknown(E, ell)
    try
        return admitsRationalIsogeny(E, ell);
    catch err
        return true;
    end try;
end function;

// Phase 2, shared by the prime and prime-power analytic paths: from the recognized
// rational jacobian invariants recIC (with source psi recPsi) and the recognized
// rational product j-pairs, reconstruct and twist-pin Q-models, canonicalize, and
// match a source psi to each returned curve. A twist-validated quotient is emitted; a
// bielliptic / Mestre one is a genuine Galois-stable gluing counted but not emitted; a
// look-alike (twist-rejected, order-2 automorphisms, over Q) is dropped. Returns the
// canonical curve list cs, the per-curve source psi list, and the analytic quotient
// count jacCount + #distinct products (the certificate unit).
function gluingEmitCurves(recIC, recPsi, products, E1, E2, TraceBound)
    curves := []; sourcePsis := []; jacCount := 0;
    for ic in Setseq(Seqset(recIC)) do
        okC, C := CurveFromInvariants(ic, E1, E2 : TraceBound := TraceBound);
        if okC then
            jacCount +:= 1;
            Append(~curves, C);
            for j in [1 .. #recIC] do
                if recIC[j] eq ic then Append(~sourcePsis, recPsi[j]); break; end if;
            end for;
        else
            genuine := false;
            try
                C0 := HyperellipticCurveFromIgusaClebsch(ic);
                genuine := (Type(BaseRing(C0)) ne FldRat) or (#GeometricAutomorphismGroup(C0) gt 2);
            catch e
                genuine := false;
            end try;
            if genuine then
                jacCount +:= 1;
                vprintf Gluing: "Genus2Gluings: quotient recognized but not emitted (bielliptic / Mestre), counted only\n";
            else
                vprintf Gluing: "Genus2Gluings: recognized quotient is a non-stable look-alike, dropping\n";
            end if;
        end if;
    end for;
    analyticCount := jacCount + #Seqset(products);
    cs := CanonicalGluingList(curves);
    // A gluing matrix per returned curve. CanonicalGluingList reduces, dedups by
    // isomorphism and sorts, so match each returned curve back to a source psi.
    psisOut := [];
    for c in cs do
        for i in [1 .. #curves] do
            if IsIsomorphic(c, curves[i]) then Append(~psisOut, sourcePsis[i]); break; end if;
        end for;
    end for;
    return cs, psisOut, analyticCount;
end function;

// -------- Prime-power levels n = ell^e via lifting (Task 10) --------
//
// Genus2Gluings at n = ell^e with e >= 2 does NOT enumerate every anti-symplectic psi
// mod ell^e (there are |SL2(Z/ell^e)| of them, growing like ell^(3e)). Instead it
// LIFTS: run the prime level ell, keep the psi whose quotient looks Q-rational, and
// lift to ell^2, ..., ell^e, at each step keeping only conjugation-filtered lifts whose
// quotient still looks rational at a cheap sweep precision. Target-level recognition,
// reconstruction and twist pinning then proceed exactly as at a prime.
//
// Why lifting the rational-LOOKING level-ell survivors (not just the ones that
// reconstruct to a curve) reaches every genuine level-ell^e gluing: the graph of a
// Galois-stable psi mod ell^e reduces mod ell to ell^(e-1) * graph(psi) =
// graph(psi | E1[ell]), which is Galois-stable, so its quotient is a Q-rational ppas
// and passes the near-real gate at level ell. Hence every stable level-ell^e graph sits
// above a level-ell survivor and is produced by lifting it. The converse can fail (a
// rational-looking level-ell graph need not lift to a stable one), which only costs
// extra candidates that the target-level recognition and twist pinning discard.
// Gap: quotients whose invariants are projectively I2 = 0 are not recognized by the
// I2-pivot (recognize.m's normalizeIC) and so would be missed at these uncertified
// levels (documented limitation).

// The near-real gate for the lift sweep: does this numeric quotient look defined over
// Q? Jacobian type reuses LooksRationalIC (the fixed 1e-10 I2-normalized gate of the
// prime pass 1); product type asks the two j-invariants to be near-real.
function quotientLooksRational(typ, inv)
    if typ eq "jacobian" then return LooksRationalIC(inv); end if;
    return forall{j : j in inv | Abs(Im(j)) lt 10^-10};
end function;

// One sweep level: compute the quotient of each candidate psi at modulus m against the
// period bases ws1, ws2 and keep those looking Q-rational.
function gluingSweepLevel(ws1, ws2, cand, m)
    out := [];
    for psi in cand do
        typ, inv := NumericInvariants(SmallFromBig(GluedBigPeriodMatrix(ws1, ws2, psi, m)));
        if quotientLooksRational(typ, inv) then Append(~out, psi); end if;
    end for;
    return out;
end function;

// Lift sweep from the prime level ell up to ell^e (e >= 1), at the cheap sweep precision
// of ws1, ws2. c1, c2 are the integer conjugation matrices from EllipticPeriodBasis; the
// conjugation-equivariance filter psi c1^T = c2^T psi (analytic.m) is a necessary
// condition at EVERY modulus, applied here with c1, c2 reduced mod ell^k. Returns the
// final-level survivor psi list (over Integers(ell^e)) and a flag that is false when an
// intermediate level emptied (hence no gluing at ell^e).
//
// The base level uses AntiSymplecticIsomorphisms(ell) (ModMinus default): for odd ell
// that is one representative per {psi, -psi} orbit, and since the ell^3 lifts of psi are
// disjoint from those of -psi = -(lifts of psi), lifting stays one-per-orbit up the
// tower. At ell = 2, -psi = psi mod 2 so the base carries all det -1 matrices and a
// psi's lifts to mod 4 include -psi's; the resulting {psi, -psi} redundancy at ell^e is
// harmless and collapsed later by the moduli-level Seqset dedup in gluingEmitCurves.
function gluingLiftSweep(ws1, c1, ws2, c2, ell, e)
    m := ell;
    Rm := Integers(m);
    c1T := Transpose(ChangeRing(c1, Rm)); c2T := Transpose(ChangeRing(c2, Rm));
    cand := [psi : psi in AntiSymplecticIsomorphisms(ell)
             | ChangeRing(psi, Rm) * c1T eq c2T * ChangeRing(psi, Rm)];
    survivors := gluingSweepLevel(ws1, ws2, cand, m);
    vprintf Gluing: "gluingLiftSweep: level %o^1: %o candidate(s) -> %o rational-looking survivor(s)\n",
        ell, #cand, #survivors;
    for k in [2 .. e] do
        if #survivors eq 0 then return survivors, false; end if;
        m := ell^k;
        Rm := Integers(m);
        c1T := Transpose(ChangeRing(c1, Rm)); c2T := Transpose(ChangeRing(c2, Rm));
        lifts := [];
        for psi in survivors do
            for N in LiftAntiSymplectic(psi, ell, k - 1) do
                if ChangeRing(N, Rm) * c1T eq c2T * ChangeRing(N, Rm) then Append(~lifts, N); end if;
            end for;
        end for;
        survivors := gluingSweepLevel(ws1, ws2, lifts, m);
        vprintf Gluing: "gluingLiftSweep: level %o^%o: %o conjugation-filtered lift(s) -> %o survivor(s)\n",
            ell, k, #lifts, #survivors;
    end for;
    return survivors, #survivors gt 0;
end function;

// Full-precision recognition of the final-level lift survivors (the e >= 2 analog of
// gluingAnalyticEnumerate's pass 2): recompute the period bases at the target
// precision, classify each survivor psi (modulus m = ell^e), recognize rational
// jacobian Igusa-Clebsch quadruples (recIC / recPsi) and rational product j-pairs
// (products), doubling the precision while that strictly grows the jacobian count (at
// most 3 times). Returns recIC, recPsi, products and the precision used.
function gluingRecognizeSurvivors(E1, E2, survivors, m, prec)
    doublings := 0; prevRecognized := -1;
    recIC := []; recPsi := []; products := []; usedPrec := prec;
    while true do
        ws1 := EllipticPeriodBasis(E1, prec);
        ws2 := EllipticPeriodBasis(E2, prec);
        usedPrec := fieldPrecision(Universe(ws1));
        recIC := []; recPsi := []; products := []; nearFail := 0;
        for psi in survivors do
            typ, inv := NumericInvariants(SmallFromBig(GluedBigPeriodMatrix(ws1, ws2, psi, m)));
            if typ eq "product" then
                ok1, j1 := RecognizeRational(inv[1]);
                ok2, j2 := RecognizeRational(inv[2]);
                if ok1 and ok2 then Append(~products, <j1, j2>); end if;
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
                vprintf Gluing: "Genus2Gluings: %o near-rational lift quotient(s) unrecognized, dropping\n", nearFail;
            end if;
            break;
        end if;
        prevRecognized := #recIC; prec := 2 * usedPrec; doublings +:= 1;
    end while;
    return recIC, recPsi, products, usedPrec;
end function;

// Level-n recognition for the composite path (gluingCompositeCRT). Same shape as
// gluingRecognizeSurvivors but the precision-doubling watches BOTH the jacobian and the
// product recognition counts: a composite (n, n)-gluing is often product-type (a stable graph
// whose quotient decomposes, e.g. 54a1 x 54b1 at n = 6), and its two near-real j-invariants can
// need more digits than the level-n heuristic gives, so a jacobian-only retry (as in
// gluingRecognizeSurvivors, which the e >= 2 prime-power path relies on unchanged) would break
// out immediately with the products undercounted. Doubles while the total recognized count
// (jacobian + product tuples, no {psi, -psi} dedup: this is the GRAPH count the certificate
// compares) strictly grows and a near-real quotient stays unrecognized, at most 3 times.
function gluingRecognizeComposite(E1, E2, survivors, m, prec)
    doublings := 0; prevRec := -1;
    recIC := []; recPsi := []; products := []; usedPrec := prec;
    while true do
        ws1 := EllipticPeriodBasis(E1, prec);
        ws2 := EllipticPeriodBasis(E2, prec);
        usedPrec := fieldPrecision(Universe(ws1));
        gate := 10^(-(usedPrec div 2));
        recIC := []; recPsi := []; products := []; nearFail := 0;
        for psi in survivors do
            typ, inv := NumericInvariants(SmallFromBig(GluedBigPeriodMatrix(ws1, ws2, psi, m)));
            if typ eq "product" then
                ok1, j1 := RecognizeRational(inv[1]);
                ok2, j2 := RecognizeRational(inv[2]);
                if ok1 and ok2 then
                    Append(~products, <j1, j2>);
                elif forall{j : j in inv | Abs(Im(j)) lt gate} then
                    nearFail +:= 1;   // near-real product j's unrecognized: a precision shortfall
                end if;
                continue;
            end if;
            okIC, ICQ := RecognizeIgusaClebsch(inv);
            if okIC then
                Append(~recIC, ICQ); Append(~recPsi, psi);
            elif IgusaClebschNearRational(inv) then
                nearFail +:= 1;
            end if;
        end for;
        totalRec := #recIC + #products;
        if nearFail eq 0 or doublings ge 3 or totalRec le prevRec then
            if nearFail gt 0 then
                vprintf Gluing: "gluingCompositeCRT: %o near-real quotient(s) unrecognized, dropping\n", nearFail;
            end if;
            break;
        end if;
        prevRec := totalRec; prec := 2 * usedPrec; doublings +:= 1;
    end while;
    return recIC, recPsi, products, usedPrec;
end function;

// -------- Composite levels n = prod ell_i^e_i via CRT composition (Task 11) --------
//
// (Z/n)^2 = direct sum_i (Z/ell_i^e_i)^2 by CRT, so an anti-symplectic psi mod n is exactly
// a CRT tuple of anti-symplectic psi_i mod ell_i^e_i (det psi = -1 mod n iff det psi_i = -1
// mod ell_i^e_i for every i). E[n] = direct sum_i E[ell_i^e_i] as Galois modules, so
// graph(psi) is Galois-stable iff every block graph(psi_i) is stable. Hence: sweep each
// prime-power block for its rational-looking survivor psi_i (gluingLiftSweep, the same
// per-block pipeline the prime-power path uses), CRT every combination of block survivors
// (CRTAntiSymplectic) into a psi mod n, run the analytic quotient at level n on those
// combinations only, and finish through the recognition/reconstruction/twist phase
// (gluingRecognizeComposite + the shared gluingEmitCurves). The level-n recognition is the precise gate:
// a CRT tuple that mixes any block's rational-LOOKING-but-not-stable look-alike is not
// Galois-stable at n, so its level-n quotient is not Q-rational and drops out, leaving exactly
// the genuine stable gluings (per-block sweeps are kept full both-sign so no combination is
// missed).
//
// Certificate, at the GRAPH level. The exact layer (exact.m) runs at each block's PRIME
// level. A composite graph is stable iff every block graph is, so the composite stable GRAPH
// count is the PRODUCT of the per-block stable graph counts. Graph, not quotient: the
// {psi, -psi} identification (-1 mod n negates every block at once) means orbit counts do NOT
// multiply across several odd blocks (|tuples| = prod graphs_i, orbits = |tuples|/2 for n > 2,
// so prod(orbits_i) undercounts by 2^(#oddblocks - 1)); multiplying at the graph level
// sidesteps that. The analytic side counts graphs too: one recognized rational tuple = one
// graph, no {psi, -psi} dedup (#recIC jacobian + #products product). certified iff every block
// certified AND the analytic graph count equals the product of block graph counts; a
// disagreement is the same hard error as the prime certificate, ell replaced by n. If any
// block is certified empty (congruence obstruction on its ell^e-part, exact count 0, or the
// e >= 2 certified-empty-by-reduction) the composite is certified empty and returns early.
// Composite products are NOT rational-isogeny gated (unlike the prime path): the level-n
// recognition already discards non-stable graphs, and a genuine product block-gluing need not
// give either curve a rational ell-isogeny; a residual look-alike would surface as a loud
// certificate mismatch, never a silent false certification.
function gluingCompositeCRT(E1, E2, n, PrecParam, Proof, TraceBound)
    fact := Factorization(n);
    k := #fact;
    ells := [f[1] : f in fact];
    es := [f[2] : f in fact];
    moduli := [ells[i]^es[i] : i in [1 .. k]];
    strict := (Proof cmpeq true);
    // The exact certificate assumes Aut(E1 x E2) = {+-1}^2. A Q-isogeny (Hom over Q) enlarges
    // it, as the prime path already documents; so does a GEOMETRIC one, which for the block
    // counts to compose must be excluded here too. A shared j-invariant is the common witness
    // (E1, E2 isomorphic over Qbar, the twist-family degeneracy, e.g. 54a1 x 54b1, j = 9261/8):
    // the CRT of stable BLOCK graphs stays valid, but such a graph's level-n quotient need not
    // be a recognizable rational curve, so prod(block graph counts) over-predicts the analytic
    // count. These pairs are left traces-only rather than risking a spurious mismatch abort.
    degenerate := IsIsogenous(E1, E2) or (jInvariant(E1) eq jInvariant(E2));

    // Congruence prefilter first: a rational n-gluing forces n | GluingModulus, and n | N iff
    // every ell_i^e_i | N. When N is conclusive and n does not divide it, each block whose
    // ell^e-part fails to divide N is congruence-certified empty, so the composite is.
    N, inconclusive := GluingModulus(E1, E2);
    if (not inconclusive) and (N mod n ne 0) then
        blocks := [];
        for i in [1 .. k] do
            if N mod moduli[i] ne 0 then
                Append(~blocks, <ells[i], es[i], 0, 0, true>);
            else
                Append(~blocks, <ells[i], es[i], -1, -1, false>);
            end if;
        end for;
        info := rec< gluingInfoFmt() | n := n, proof := "certified",
            blocks := blocks, psis := [], products := [],
            precision := 0, tracebound := TraceBound >;
        return emptyCurves(), info;
    end if;

    // Per-block: exact certified-empty short-circuit + rational-looking survivor sweep.
    blockTuples := [];
    blockGraphCounts := [];
    blockSurvivorsFull := [* *];
    allCertified := true;
    for i in [1 .. k] do
        ell := ells[i]; e := es[i];
        // Mirror the prime path: exact certifies "Auto" blocks up to ell = 13, true every
        // prime, and never runs for isogenous pairs (Aut(E1 x E2) then exceeds {+-1}^2).
        runExactBlock := ((Proof cmpeq true) or (Proof cmpeq "Auto" and ell le 13)) and not degenerate;
        stableOrbit := -1; graphCount := -1; certified := false;
        if runExactBlock then
            so, srec := GaloisStableGluings(E1, E2, ell);
            if so lt 0 then
                error if strict,
                    Sprintf("gluing certificate unavailable at ell=%o: exact layer declined (Galois group order %o over the degree bound)", ell, srec`group_order);
            elif so eq 0 then
                // Certified empty (e = 1 directly; e >= 2 by reduction mod ell). The whole
                // composite is certified empty: return with this block as the witness.
                blocks := blockTuples cat [<ell, e, 0, 0, true>];
                for j in [i + 1 .. k] do Append(~blocks, <ells[j], es[j], -1, -1, false>); end for;
                info := rec< gluingInfoFmt() | n := n, proof := "certified",
                    blocks := blocks, psis := [], products := [],
                    precision := 0, tracebound := TraceBound >;
                return emptyCurves(), info;
            elif e eq 1 then
                stableOrbit := so; graphCount := srec`graph_count; certified := true;
            end if;   // e >= 2 with so > 0: no exact ell^e layer, block stays uncertified
        end if;
        if not certified then allCertified := false; end if;

        sweepPrec := GluingPrecisionHeuristic(E1, E2, ell);
        ws1s, cc1 := EllipticPeriodBasis(E1, sweepPrec);
        ws2s, cc2 := EllipticPeriodBasis(E2, sweepPrec);
        survivors := gluingLiftSweep(ws1s, cc1, ws2s, cc2, ell, e);
        // Full both-sign survivor list (graph level). gluingLiftSweep returns one per
        // {psi, -psi} orbit; adjoin negatives so every combination of block graphs is CRT'd
        // (disjoint for a sign-nontrivial block, ell odd or e >= 2; identical mod 2 at e = 1).
        full := [* *]; seen := {};
        for psi in survivors do
            for M in [psi, -psi] do
                key := Eltseq(M);
                if key notin seen then Include(~seen, key); Append(~full, M); end if;
            end for;
        end for;
        Append(~blockSurvivorsFull, full);
        Append(~blockGraphCounts, graphCount);
        Append(~blockTuples, <ell, e, stableOrbit, #survivors, certified>);
        vprintf Gluing: "gluingCompositeCRT: block %o^%o: %o survivor orbit(s) -> %o graph(s); exact graph %o, certified %o\n",
            ell, e, #survivors, #full, graphCount, certified;
    end for;

    // CRT-combine every tuple of block survivors into a psi mod n.
    curPsis := blockSurvivorsFull[1]; curMod := moduli[1];
    for i in [2 .. k] do
        newPsis := [* *];
        for a in curPsis do
            for b in blockSurvivorsFull[i] do
                Append(~newPsis, CRTAntiSymplectic([* a, b *], [curMod, moduli[i]]));
            end for;
        end for;
        curPsis := newPsis; curMod := curMod * moduli[i];
    end for;
    vprintf Gluing: "gluingCompositeCRT: %o CRT combination(s) at level %o\n", #curPsis, n;

    // Level-n recognition (product-aware precision-doubling recognizer) + reconstruction.
    prec := (PrecParam cmpeq false) select GluingPrecisionHeuristic(E1, E2, n) else PrecParam;
    recIC, recPsi, products, usedPrec := gluingRecognizeComposite(E1, E2, curPsis, n, prec);
    analyticGraph := #recIC + #products;
    cs, psisOut, _ := gluingEmitCurves(recIC, recPsi, products, E1, E2, TraceBound);

    proof := "traces-only";
    if allCertified then
        exactGraph := &*blockGraphCounts;   // all blocks certified => every graphCount >= 0
        error if analyticGraph ne exactGraph,
            Sprintf("gluing certificate mismatch at n=%o: exact %o vs analytic %o", n, exactGraph, analyticGraph);
        proof := "certified";
    end if;
    info := rec< gluingInfoFmt() | n := n, proof := proof,
        blocks := blockTuples, psis := psisOut, products := products,
        precision := usedPrec, tracebound := TraceBound >;
    return cs, info;
end function;

intrinsic Genus2Gluings(E1::CrvEll, E2::CrvEll, n::RngIntElt
    : Algorithm := "Auto", Precision := false, Proof := "Auto", TraceBound := 1000)
    -> SeqEnum, Rec
{The genus-2 curves over Q gluing E1 and E2 along full n-torsion, for any n >= 2, with
a GluingInfo record. A composite n factors into prime-power blocks ell_i^e_i: each block
is swept for its rational-looking survivor psi_i, every CRT combination psi = (psi_i) mod
n is realized analytically at level n, and the shared recognition/twist phase reconstructs
the curves (see the gluingCompositeCRT header). Algorithm is one of "Auto" (BHLS closed
formulas when n is 2 or 3 and the curves are non-isomorphic, else analytic periods),
"Algebraic" (BHLS, n = 2 or 3 only), or "Periods" (always analytic). For e >= 2 the
analytic path lifts a level-ell survivor set up the tower (see the gluingLiftSweep header)
rather than enumerating psi mod ell^e. Precision sets the STARTING analytic precision
(default the analytic precision heuristic); a near-rational recognition failure retries
at double the precision (see the Pipeline note above). TraceBound bounds the Euler-factor
twist certificate. Proof drives the exact
completeness certificate (exact.m, GaloisStableGluings): "Auto" certifies prime blocks up
to ell = 13, true certifies every prime and errors if the exact layer declines, false
skips it. The exact layer runs at each block's PRIME level only; for e >= 2 blocks it is
out of scope in v1 (division polynomials at ell^e are too large), so those blocks are
"traces-only" EXCEPT that a prime level certified empty forces ell^e empty by reduction
(certified). The same shortcut applies to a composite n (gluingCompositeCRT): one
certified-empty block (congruence-obstructed or exact-certified-empty) proves the whole
composite empty, and the composite returns "certified" immediately with the remaining
blocks in the tuple reported as unexamined sentinels, not independently certified. The
GluingInfo fields are n, proof ("certified" when every block certified,
else "traces-only"), blocks (one <ell, e, stable_count, analytic_count, certified> tuple
per prime-power block; stable_count is the exact Galois-stable QUOTIENT count when
certified, -1 otherwise (including every non-empty e >= 2 block), analytic_count the
block's rational-looking survivor orbit count; the composite certificate compares the
PRODUCT of block graph counts against the analytic graph count and a disagreement is a hard
error), psis (a gluing matrix per returned curve), products (recognized <j1, j2> pairs of
the product-type quotients), precision (digits of the successful analytic pass), and
tracebound. Base field is Q, or (experimental, Task 13) a number field K with E1, E2 over the
same K and the first infinite place of K real: number-field inputs run the analytic periods
path at one fixed real embedding (EllipticPeriodBasis), recognize over Q, and are always
"traces-only" (the exact/BHLS layers are Q-only) except for congruence-certified-empty
blocks; prime-power levels only.}
    require n ge 2: "n must be at least 2";
    require Algorithm in ["Auto", "Algebraic", "Periods"]:
        "Algorithm must be one of \"Auto\", \"Algebraic\", \"Periods\"";
    K := BaseRing(E1);
    require BaseRing(E2) cmpeq K: "E1 and E2 must share a base field";
    require Type(K) eq FldRat or ISA(Type(K), FldNum):
        "Genus2Gluings requires E1, E2 over Q or a common number field";
    isNF := Type(K) ne FldRat;
    require Proof cmpeq "Auto" or Proof cmpeq true or Proof cmpeq false:
        "Proof must be \"Auto\", true, or false";
    require not (isNF and Proof cmpeq true):
        "Proof := true is unavailable over number fields; the exact layer is Q-only in v1";

    // Composite n = prod ell_i^e_i: CRT composition of prime-power blocks (gluingCompositeCRT
    // header). "Algebraic" is undefined off n in {2, 3}, so reject it here (matching the
    // prime-power require below) before the analytic composition runs.
    if not IsPrimePower(n) then
        require Algorithm ne "Algebraic":
            "the \"Algebraic\" algorithm is only defined for n in {2, 3}";
        require not isNF:
            "number-field inputs (Task 13, experimental) support prime-power levels only; composite n is not yet handled over a number field";
        return gluingCompositeCRT(E1, E2, n, Precision, Proof, TraceBound);
    end if;

    _, ell, e := IsPrimePower(n);

    // Exact completeness certificate (exact.m). "Auto" certifies prime blocks up
    // to ell = 13 (the exact layer's practical reach at DegreeBound 2400); true
    // certifies every prime and errors if the exact layer declines; false skips.
    // The certificate is only run for non-isogenous E1, E2: when Hom(E1, E2) != 0
    // the surface E1 x E2 carries extra endomorphisms, so Aut(E1 x E2) is larger
    // than {+-1}^2 and product-type quotients proliferate (a cross-isogeny already
    // gives (E1 x E2)/graph(psi) = E1 x E2). The analytic side then cannot
    // reconcile its rational-quotient count with the exact graph count by the
    // {psi, -psi}-fold and rational-isogeny gate this module uses, so isogenous
    // pairs (the degenerate iso/eq/cm inputs) stay "traces-only".
    // Over a number field the exact completeness layer (exact.m) and BHLS both
    // require Q, and IsIsogenous is not available for CrvEll[FldNum], so force the
    // certificate off: number-field blocks are traces-only in v1 (the congruence
    // certified-empty path below still applies, GluingModulus holding over K). The
    // "and" short-circuits, so IsIsogenous is never evaluated over K.
    runExact := (not isNF)
        and ((Proof cmpeq true) or (Proof cmpeq "Auto" and ell le 13))
        and not IsIsogenous(E1, E2);
    strict := (not isNF) and (Proof cmpeq true);

    // Dispatch validation, ahead of the fast path below: an Algorithm choice
    // that does not apply to this (E1, E2, n) must error, even when the pair
    // is certifiably empty (e.g. Algorithm := "Algebraic" at n = 5 must raise
    // the n in {2, 3} require, not silently return an empty list).
    tryAlgebraic := false;
    algebraicStrict := false;   // "Algebraic": let the BHLS requires propagate
    if Algorithm eq "Algebraic" then
        require not isNF:
            "the \"Algebraic\" algorithm (BHLS) requires E1, E2 over Q; number-field inputs use the analytic periods path";
        require n in {2, 3}: "the \"Algebraic\" algorithm is only defined for n in {2, 3}";
        require not IsIsomorphic(E1, E2):
            "the \"Algebraic\" algorithm requires non-isomorphic curves";
        tryAlgebraic := true; algebraicStrict := true;
    elif Algorithm eq "Auto" then
        // BHLS is Q-only; over K the "and" short-circuits before IsIsomorphic.
        if (not isNF) and n in {2, 3} and not IsIsomorphic(E1, E2) then tryAlgebraic := true; end if;
    end if;

    // Congruence-obstruction fast path: a rational n-gluing forces n | a_p(E1) -
    // a_p(E2) at every good p, hence n | GluingModulus. When that modulus is
    // conclusive (not every scanned trace agreed) and n does not divide it, the
    // gluing set is provably empty with no analytic work.
    N, inconclusive := GluingModulus(E1, E2);
    if (not inconclusive) and (N mod n ne 0) then
        // Congruence obstruction: a rigorous, independent proof of emptiness (n
        // does not divide a good-prime trace difference), so 0 stable graphs and
        // 0 rational quotients, block certified without the exact layer. Necessary
        // at prime powers too: a rational ell^e-gluing forces a_p(E1) = a_p(E2)
        // mod ell^e, hence ell^e | GluingModulus.
        info := rec< gluingInfoFmt() | n := n, proof := "certified",
            blocks := [<ell, e, 0, 0, true>], psis := [], products := [],
            precision := 0, tracebound := TraceBound >;
        return emptyCurves(), info;
    end if;

    // Prime-power levels n = ell^e with e >= 2: lift a level-ell survivor set rather
    // than enumerate every psi mod ell^e (see the gluingLiftSweep header). No algebraic
    // formulas apply (BHLS is n in {2, 3}, hence e = 1), so this is the only path, and
    // the exact completeness layer is out of scope here (v1), so blocks are traces-only
    // save the certified-empty upgrade below.
    if e ge 2 then
        // Certified-empty-by-reduction: a Galois-stable ell^e graph reduces mod ell to a
        // Galois-stable ell graph, so if the exact layer certifies the prime level empty,
        // ell^e is provably empty for every e. runExact already excludes isogenous pairs
        // and respects Proof; GaloisStableGluings declines (returns -1, not 0) when the
        // Galois group is too big, in which case we fall through to the sweep.
        if runExact and GaloisStableGluings(E1, E2, ell) eq 0 then
            info := rec< gluingInfoFmt() | n := n, proof := "certified",
                blocks := [<ell, e, 0, 0, true>], psis := [], products := [],
                precision := 0, tracebound := TraceBound >;
            return emptyCurves(), info;
        end if;

        // Lift sweep at the cheap prime-level sweep precision; recognition then reruns on
        // the survivors at full target precision.
        sweepPrec := GluingPrecisionHeuristic(E1, E2, ell);
        ws1s, cc1 := EllipticPeriodBasis(E1, sweepPrec);
        ws2s, cc2 := EllipticPeriodBasis(E2, sweepPrec);
        survivors, nonempty := gluingLiftSweep(ws1s, cc1, ws2s, cc2, ell, e);
        if not nonempty then
            // An intermediate level had no rational-looking survivor, so no ell^e gluing.
            // Not certified (the sweep is numeric; the prime level was not certified empty
            // above), so traces-only.
            info := rec< gluingInfoFmt() | n := n, proof := "traces-only",
                blocks := [<ell, e, -1, 0, false>], psis := [], products := [],
                precision := sweepPrec, tracebound := TraceBound >;
            return emptyCurves(), info;
        end if;
        prec := (Precision cmpeq false) select GluingPrecisionHeuristic(E1, E2, n) else Precision;
        recIC, recPsi, products, usedPrec := gluingRecognizeSurvivors(E1, E2, survivors, n, prec);
        // ell can exceed 59 here (e >= 2 lift path), past ClassicalModularPolynomial's
        // database: admitsRationalIsogenyOrUnknown treats that as "cannot rule out".
        if not (admitsRationalIsogenyOrUnknown(E1, ell) and admitsRationalIsogenyOrUnknown(E2, ell)) then products := []; end if;
        cs, psisOut, analyticCount := gluingEmitCurves(recIC, recPsi, products, E1, E2, TraceBound);
        // No exact ell^e layer in v1, so stable_count is -1 and the block is traces-only
        // even when curves are emitted (each still passes the Euler-factor twist
        // certificate to TraceBound in CurveFromInvariants).
        info := rec< gluingInfoFmt() | n := n, proof := "traces-only",
            blocks := [<ell, e, -1, analyticCount, false>], psis := psisOut, products := products,
            precision := usedPrec, tracebound := TraceBound >;
        return cs, info;
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
        if runExact then
            // BHLS is a complete algorithm for the (n, n)-jacobian gluings at
            // n in {2, 3}, so #cs is the exact jacobian-quotient count (recognition
            // is unreliable here: at small mod-n image the conjugation filter barely
            // cuts the pool and look-alike quotients isogenous to E1 x E2 twist-
            // validate). Only the product quotients, which BHLS does not emit, are
            // recovered from the period enumeration (rational-isogeny gated to drop
            // product look-alikes; see admitsRationalIsogeny).
            prec := (Precision cmpeq false) select GluingPrecisionHeuristic(E1, E2, n) else Precision;
            _, _, products, usedPrec := gluingAnalyticEnumerate(E1, E2, n, prec);
            // n in {2, 3} here (usedAlgebraic implies BHLS domain), well inside
            // ClassicalModularPolynomial's range; only worth the gate when there is a
            // product quotient to gate.
            if #products gt 0 and not (admitsRationalIsogeny(E1, n) and admitsRationalIsogeny(E2, n)) then products := []; end if;
            analyticCount := #cs + #Seqset(products);
            block, blockProof := GluingCertificateBlock(E1, E2, n, analyticCount, strict);
        else
            products := []; usedPrec := 0;
            block := <n, 1, -1, #cs, false>; blockProof := "traces-only";
        end if;
        info := rec< gluingInfoFmt() | n := n, proof := blockProof,
            blocks := [block], psis := [], products := products,
            precision := usedPrec, tracebound := TraceBound >;
        return cs, info;
    end if;

    // Analytic (periods) path. The two-pass precision sweep enumerating and
    // recognizing the anti-symplectic quotients is gluingAnalyticEnumerate above;
    // here we reconstruct Q-models from the recognized jacobian invariants and run
    // the completeness certificate.
    prec := (Precision cmpeq false) select GluingPrecisionHeuristic(E1, E2, n) else Precision;
    recIC, recPsi, products, usedPrec := gluingAnalyticEnumerate(E1, E2, n, prec);
    // n is an arbitrary prime here (the general periods path) and can exceed 59, past
    // ClassicalModularPolynomial's database: admitsRationalIsogenyOrUnknown treats that
    // as "cannot rule out".
    if not (admitsRationalIsogenyOrUnknown(E1, n) and admitsRationalIsogenyOrUnknown(E2, n)) then products := []; end if;

    // Phase 2 (reconstruction + count), shared with the prime-power path via
    // gluingEmitCurves: one pass over the distinct recognized jacobian moduli. A
    // twist-validated quotient is emitted as a curve; a bielliptic / Mestre one is a
    // genuine Galois-stable gluing counted but not emitted; a look-alike (twist-rejected,
    // order-2 automorphisms, over Q) has no rational isogeny to E1 x E2 and is dropped
    // from both. analyticCount is the genuine jacobian quotients plus the
    // rational-isogeny-gated products, the same unit as the exact layer.
    cs, psisOut, analyticCount := gluingEmitCurves(recIC, recPsi, products, E1, E2, TraceBound);

    // The analytic recognition is a reliable rational-quotient count only at odd n,
    // where the conjugation filter is tight; at n = 2 it barely cuts the candidate
    // pool and look-alike jacobians isogenous to E1 x E2 twist-validate, so the n = 2
    // certificate is left to the BHLS path (which Auto always takes at n = 2, non-
    // isomorphic). Under Algorithm := "Periods" an n = 2 block stays "traces-only".
    if runExact and n ge 3 then
        block, blockProof := GluingCertificateBlock(E1, E2, n, analyticCount, strict);
    else
        block := <n, 1, -1, analyticCount, false>; blockProof := "traces-only";
    end if;
    info := rec< gluingInfoFmt() | n := n, proof := blockProof,
        blocks := [block], psis := psisOut, products := products,
        precision := usedPrec, tracebound := TraceBound >;
    return cs, info;
end intrinsic;

intrinsic Genus2Gluings(f1::RngUPolElt, f2::RngUPolElt, n::RngIntElt
    : Algorithm := "Auto", Precision := false, Proof := "Auto", TraceBound := 1000)
    -> SeqEnum, Rec
{Genus2Gluings for the elliptic curves y^2 = f1, y^2 = f2 (f1, f2 monic univariate
cubics over Q).}
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
