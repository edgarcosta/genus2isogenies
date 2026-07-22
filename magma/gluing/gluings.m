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
 * The exact completeness certificate (exact.m) closes the loop, on the paper's
 * Algorithm 3.4 scope only: geometrically NONisogenous pairs, Hom_Qbar(E1, E2) = 0
 * (the clean locus). There GaloisStableGluings counts the Galois-stable
 * anti-symplectic gluings exactly (as M/-M-orbit quotients) and
 * GluingCertificateBlock compares that against the analytic count of distinct
 * emitted rational jacobian quotients (#Seqset over recIC, each reconstructed AND
 * twist-pinned by CurveFromInvariants). On the clean locus those units agree:
 * Kani's reducibility criterion makes a product-type graph quotient impossible (it
 * would force a nonzero geometric isogeny E1 -> E2 of degree k(n-k)), so a
 * recognized product is always a look-alike and never counts (reported in
 * info`products only); and a recognized quotient rejected by CurveFromInvariants
 * never counts either, it leaves the block UNCERTIFIABLE (traces-only), because
 * rational moduli plus extra automorphisms do not imply a Galois-stable graph (a
 * Galois-stable automorphism ORBIT of graphs already gives rational moduli with no
 * stable member). Under Proof := "Auto" (prime n <= 13) or true, a count match
 * certifies the block and a disagreement is a hard error; the proof field is
 * "certified" iff every block certified. A geometrically isogenous pair
 * (GeometricallyIsogenous, modulus.m; the paper's still-in-development "isogenous
 * to a square" case) is NEVER certified: proof stays "traces-only" and, at prime
 * blocks where the certificate would have run, the exact stable count is still
 * computed and recorded in the block tuple as diagnostic data (certified false,
 * never compared, never an error).
 */

function gluingInfoFmt()
    return recformat< n : RngIntElt, proof : MonStgElt, blocks : SeqEnum,
                      psis : SeqEnum, products : SeqEnum,
                      precision : RngIntElt, tracebound : RngIntElt >;
end function;

function emptyCurves()
    return [Parent(HyperellipticCurve(PolynomialRing(Rationals()).1^5 + 1)) | ];
end function;

// Strict Proof := true (option (a)) belt-and-braces: no successful return may carry a
// traces-only certificate. The site-specific errors at each decision point are the primary
// gate (they give an actionable message and fire before the return); this is the safety net
// that catches any traces-only return left unguarded by a future edit. It checks info`proof,
// NOT the per-block certified flags: a certified-empty return legitimately carries uncertified
// sentinel blocks (the composite short-circuit records later blocks as <ell, e, -1, -1, false>)
// under an overall proof of "certified", which a block-level check would wrongly reject.
procedure strictGuard(strict, info)
    error if strict and info`proof eq "traces-only",
        "Proof := true requires a non-traces-only certificate for every block, but this result is traces-only; use Proof := \"Auto\"";
end procedure;

// Working precision of a complex field. A file-local helper so it reaches the
// builtin Precision, which the Genus2Gluings "Precision" parameter shadows
// inside that intrinsic's body.
function fieldPrecision(C)
    return Precision(C);
end function;

// Pass 1 of the two-pass sweep (Genus2Gluings below): a coarse, fixed-
// threshold look at whether a Jacobian-type quotient is plausibly Q-rational,
// cheap enough to run on every conjugation-filtered candidate at the fixed
// low pass-1 precision (40 digits). Mirrors recognize.m's multichart dispatcher
// (the generic I2 pivot when I2 is projectively nonzero, else the I2 = 0
// weight-0 ratio charts) but with a fixed gate (1e-10) rather than
// RecognizeIgusaClebsch/IgusaClebschNearRational's precision-scaled one (1e-20 at
// 40 digits): pass 1's job is to not miss a genuine rational quotient to roundoff
// in a deliberately cheap computation, not to make the final call, so it uses the
// looser threshold. The chart selectors here (I2^5/I10, I4^5/I10^2, I6^5/I10^3)
// are themselves weight-0, hence already scale-free without the projective scale S.
function LooksRationalIC(IC)
    I2 := IC[1]; I4 := IC[2]; I6 := IC[3]; I10 := IC[4];
    gate := 10^-10;
    if Abs(I10) le gate then return false; end if;
    if Abs(I2^5 / I10) ge gate then
        return Abs(Im(I4 / I2^2)) lt gate and Abs(Im(I6 / I2^3)) lt gate and Abs(Im(I10 / I2^5)) lt gate;
    end if;
    // I2 projectively zero: the weight-0 ratio charts (I4 chart t1, t2; then the
    // I4 = 0 chart u; then the all-zero point [0,0,0,1], trivially rational).
    if Abs(I4^5 / I10^2) ge gate then
        return Abs(Im(I4 * I6 / I10)) lt gate and Abs(Im(I4^5 / I10^2)) lt gate;
    end if;
    if Abs(I6^5 / I10^3) ge gate then
        return Abs(Im(I6^5 / I10^3)) lt gate;
    end if;
    return true;
end function;

// The shared near-real gate for a numeric quotient looking defined over Q, used by the
// prime pass 1 (gluingAnalyticEnumerate) and the lift sweep (gluingSweepLevel). Jacobian
// type reuses LooksRationalIC (the fixed 1e-10 gate of the prime pass 1); product type
// asks the two j-invariants to be near-real. A "product-unsplit" quotient (theta-
// vanishing but non-diagonal reduced tau, empty invariants) is UNCLASSIFIED, so it is
// kept as a survivor: the working-precision recognizer must see it, retry it at higher
// precision, and leave the block uncertifiable if it persists; dropping it here would
// fake certainty about a quotient nothing classified.
function quotientLooksRational(typ, inv)
    if typ eq "jacobian" then return LooksRationalIC(inv); end if;
    if typ eq "product-unsplit" then return true; end if;
    return forall{j : j in inv | Abs(Im(j)) lt 10^-10};
end function;

// The single recognizer for every analytic path (prime pass 2, prime-power lift sweep, and
// composite CRT): recompute the period bases at the working precision, classify each survivor
// psi at modulus m, recognize rational jacobian Igusa-Clebsch quadruples (recIC / recPsi) and
// rational product j-pairs (products), and return them with the precision used. The precision-
// doubling watches BOTH the jacobian and the product recognition counts, since a composite
// (n, n)-gluing is often product-type (a stable graph whose quotient decomposes, e.g. 54a1 x
// 54b1 at n = 6) and its two near-real j-invariants can need more digits than the level-n
// heuristic gives; nearFail counts a near-real quotient of either kind that stayed
// unrecognized, plus every "product-unsplit" quotient (theta-vanishing, non-diagonal
// reduced tau: unclassified, so a retry at more digits is the only move available).
// It doubles while doubling HELPS, i.e. while a doubling strictly shrinks nearFail: the
// corpus-typical population of near-real yet irrational quotients (e.g. conjugate products
// over a real quadratic field, whose j's sit under the near-real gate forever) keeps
// nearFail constant, and unconditional exhaustion would then pay the 2x/4x/8x theta cost
// on essentially every prime-level call (measured ~8x overall on the corpus, up to ~24x on
// some entry classes). So the loop stops when a doubling leaves nearFail without strict
// decrease, or the budget (3 doublings) is spent. Known recall gap: a genuinely rational
// quotient needing MORE than one doubling of headroom, with nothing else recognized at the
// intermediate step, is dropped where exhaustion would have reached it. On certified paths
// such a drop surfaces as a loud count mismatch; traces-only paths (e >= 2 prime powers,
// uncertified composites) carry no such net, so a suspected drop there is addressed by
// passing a larger Precision. recIC carries no {psi, -psi} dedup: #recIC is the jacobian
// GRAPH count the composite certificate compares (products never count, see the file
// header). The fifth return value is the number of quotients STILL "product-unsplit"
// after the retry budget: each is an unclassified quotient, so a caller may not certify
// a block whose recognition left it nonzero.
function gluingRecognizeAtLevel(E1, E2, survivors, m, prec)
    doublings := 0; prevNearFail := -1;
    recIC := []; recPsi := []; products := []; usedPrec := prec; unsplit := 0;
    while true do
        ws1 := EllipticPeriodBasis(E1, prec);
        ws2 := EllipticPeriodBasis(E2, prec);
        usedPrec := fieldPrecision(Universe(ws1));
        gate := 10^(-(usedPrec div 2));
        recIC := []; recPsi := []; products := []; nearFail := 0; unsplit := 0;
        for psi in survivors do
            typ, inv := NumericInvariants(SmallFromBig(GluedBigPeriodMatrix(ws1, ws2, psi, m)));
            if typ eq "product-unsplit" then
                unsplit +:= 1; nearFail +:= 1;   // unclassified: more digits is the only move
                continue;
            end if;
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
        if nearFail eq 0 or doublings ge 3 or (prevNearFail ge 0 and nearFail ge prevNearFail) then
            if nearFail gt 0 then
                vprintf Gluing: "gluingRecognizeAtLevel: %o near-real quotient(s) unrecognized after %o doubling(s), dropping\n", nearFail, doublings;
            end if;
            if unsplit gt 0 then
                vprintf Gluing: "gluingRecognizeAtLevel: %o product-type quotient(s) still unsplit (non-diagonal reduced tau) at %o digits: unclassified, so this block cannot certify\n", unsplit, usedPrec;
            end if;
            break;
        end if;
        prevNearFail := nearFail; prec := 2 * usedPrec; doublings +:= 1;
    end while;
    return recIC, recPsi, products, usedPrec, unsplit;
end function;

// The two-pass precision sweep for the prime analytic dispatch (and the completeness
// certificate, which needs only the count). Pass 1, at a fixed cheap 40 digits, keeps every
// conjugation-filtered quotient GluedPeriodMatrices returns that looks Q-rational
// (quotientLooksRational: jacobian survivors by the I2-normalized gate, product survivors by the
// near-real j gate) and forwards their psi; pass 2 is the shared gluingRecognizeAtLevel, which
// recomputes those at working precision and recognizes the rational jacobian Igusa-Clebsch list
// recIC (with source psi recPsi) and rational product j-pairs products, doubling to exhaust the
// retry budget. Products are therefore recognized at working precision like the composite path,
// not once at the fixed pass-1 precision. Returns recIC, recPsi, products, the precision used,
// and the residual product-unsplit count (gluingRecognizeAtLevel).
//
// The certificate's analytic count is built downstream by gluingEmitCurves: DISTINCT
// recognized-and-emitted jacobian moduli (Seqset dedup, so psi and -psi count once), the same
// unit as the exact side's M/-M-orbit count. Products never count (clean-locus Kani vacuity,
// see the file header), and a recognized quotient that fails reconstruction or twist pinning
// makes the block uncertifiable instead of being counted.
function gluingAnalyticEnumerate(E1, E2, n, prec)
    pass1 := GluedPeriodMatrices(E1, E2, n : Precision := 40);
    survivors := [];
    njac := 0; nprod := 0;
    for r in pass1 do
        if r`type eq "jacobian" then njac +:= 1; else nprod +:= 1; end if;
        if quotientLooksRational(r`type, r`invariants) then Append(~survivors, r`psi); end if;
    end for;
    vprintf Gluing: "Genus2Gluings: pass 1 (precision 40) kept %o/%o candidate(s) (%o jacobian-type, %o product-type seen)\n",
        #survivors, njac + nprod, njac, nprod;
    if #survivors eq 0 then return [], [], [], 40, 0; end if;
    return gluingRecognizeAtLevel(E1, E2, survivors, n, prec);
end function;

// Phase 2, shared by the prime and prime-power analytic paths: from the recognized
// rational jacobian invariants recIC (with source psi recPsi), reconstruct and
// twist-pin Q-models, canonicalize, and match a source psi to each returned curve.
// A twist-validated quotient is emitted and counted; a quotient CurveFromInvariants
// rejects (bielliptic / Mestre reconstruction obstruction, or a twist-rejected
// look-alike) is NOT counted and instead raises the dropped flag, which makes the
// enclosing block uncertifiable. There is deliberately no "counted but not emitted"
// middle ground: rational moduli plus a geometric automorphism group of order > 2 do
// NOT imply a Galois-stable graph (a Galois-stable automorphism ORBIT of graphs
// already forces rational moduli with no stable member, so that inference
// over-counts), and this path cannot match a rejected quotient to an exact stable
// graph, so certification has to be withheld rather than guessed. Products are not
// an input: they never count (clean-locus Kani vacuity, file header). Returns the
// canonical curve list cs, the per-curve source psi list, the analytic quotient
// count (DISTINCT emitted moduli, the certificate unit), and the dropped flag.
function gluingEmitCurves(recIC, recPsi, E1, E2, TraceBound)
    curves := []; sourcePsis := []; jacCount := 0; dropped := false;
    for ic in Setseq(Seqset(recIC)) do
        okC, C := CurveFromInvariants(ic, E1, E2 : TraceBound := TraceBound);
        if okC then
            jacCount +:= 1;
            Append(~curves, C);
            for j in [1 .. #recIC] do
                if recIC[j] eq ic then Append(~sourcePsis, recPsi[j]); break; end if;
            end for;
        else
            dropped := true;
            vprintf Gluing: "Genus2Gluings: recognized quotient IC = %o rejected by CurveFromInvariants: not counted, block left uncertifiable\n", ic;
        end if;
    end for;
    analyticCount := jacCount;
    cs := CanonicalGluingList(curves);
    // A gluing matrix per returned curve. CanonicalGluingList reduces, dedups by
    // isomorphism and sorts, so match each returned curve back to a source psi.
    psisOut := [];
    for c in cs do
        for i in [1 .. #curves] do
            if IsIsomorphic(c, curves[i]) then Append(~psisOut, sourcePsis[i]); break; end if;
        end for;
    end for;
    return cs, psisOut, analyticCount, dropped;
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
// The projectively-I2 = 0 locus is now recognized by recognize.m's weight-0 ratio
// charts (icChart / recognizeICCandidates), so these lift-swept levels no longer
// drop it. As everywhere, an I2 = 0 quotient with extra automorphisms is still
// subject to the non-quadratic-twist limitation (CurveFromInvariants drops
// aut order > 2), not to a recognition gap.

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
// (gluingRecognizeAtLevel + the shared gluingEmitCurves). The level-n recognition is the precise gate:
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
// sidesteps that. The analytic side counts graphs too: one recognized rational jacobian
// tuple = one graph, no {psi, -psi} dedup (#recIC). Recognized products do NOT count: the
// certificate only runs on the clean locus (see below), where Kani's criterion rules out
// product-type graph quotients at every level n (a product would force a nonzero geometric
// isogeny E1 -> E2 of degree k(n-k)), so a recognized product is a look-alike, reported in
// info`products only. certified iff every block certified AND the analytic graph count
// equals the product of block graph counts AND the level-n recognition/reconstruction left
// nothing unclassified (no residual product-unsplit quotient, no CurveFromInvariants-
// rejected quotient); a count disagreement is the same hard error as the prime certificate,
// ell replaced by n. If any block is certified empty (congruence obstruction on its
// ell^e-part, exact count 0, or the e >= 2 certified-empty-by-reduction) the composite is
// certified empty and returns early.
function gluingCompositeCRT(E1, E2, n, PrecParam, Proof, TraceBound)
    fact := Factorization(n);
    k := #fact;
    ells := [f[1] : f in fact];
    es := [f[2] : f in fact];
    moduli := [ells[i]^es[i] : i in [1 .. k]];
    strict := (Proof cmpeq true);
    // Certificate scope (paper Algorithm 3.4; see the prime path's runExact comment): the
    // exact layer only certifies geometrically NONisogenous pairs, and IsIsogenous or
    // GeometricallyIsogenous excludes a pair outright. GeometricallyIsogenous subsumes the
    // old shared-j test (equal j is its first case) and adds the exact CM-field and
    // isogeny-class j-match cases. Computed once per call, and only when some block's
    // Proof gate could fire at all (Proof := false never needs it); excluded pairs still
    // get per-block diagnostics below.
    excluded := false;
    if (Proof cmpeq true) or ((Proof cmpeq "Auto") and Min(ells) le 13) then
        excluded := IsIsogenous(E1, E2) or GeometricallyIsogenous(E1, E2);
    end if;

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
        // prime, and never runs for excluded (Q- or geometrically isogenous) pairs.
        gateBlock := (Proof cmpeq true) or (Proof cmpeq "Auto" and ell le 13);
        runExactBlock := gateBlock and not excluded;
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
        elif gateBlock and e eq 1 then
            // Excluded pair, prime block where the certificate would have run: record the
            // exact stable count as DIAGNOSTIC data (square-case theory development). It is
            // never compared against anything and never certifies or empties a block here,
            // not even at stable count 0: the pair is outside the certificate's scope.
            stableOrbit := GaloisStableGluings(E1, E2, ell);
        end if;
        if not certified then allCertified := false; end if;

        // The sweep evaluates quotients up to the block's top level moduli[i], so its
        // precision is keyed to that top level (honoring a user Precision), not the base prime.
        sweepPrec := (PrecParam cmpeq false) select GluingPrecisionHeuristic(E1, E2, moduli[i]) else PrecParam;
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

    // Strict Proof := true (option (a)): a block left uncertified cannot enter a
    // non-traces-only composite certificate. The check is deferred to here, AFTER the loop,
    // on purpose: any block that is certified empty (congruence or exact) returns the whole
    // composite "certified" from inside the loop above, and that short-circuit must win even
    // when an EARLIER block was uncertifiable (e.g. an e >= 2 block ordered before a
    // certified-empty prime block). Reaching this point means no block was certified empty,
    // so an uncertified block (a Q-/geometrically isogenous pair, or an e >= 2 block with
    // stable gluings and no v1 exact layer) is fatal under strict.
    error if strict and not allCertified,
        Sprintf("Proof := true cannot certify n=%o: a prime-power block is uncertifiable (E1 and E2 are Q-isogenous or geometrically isogenous, or an e >= 2 block has stable gluings but no v1 exact layer above the prime level); use Proof := \"Auto\"", n);

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
    recIC, recPsi, products, usedPrec, unsplit := gluingRecognizeAtLevel(E1, E2, curPsis, n, prec);
    analyticGraph := #recIC;   // jacobian graphs only; products never count (header)
    cs, psisOut, _, dropped := gluingEmitCurves(recIC, recPsi, E1, E2, TraceBound);
    if allCertified and (unsplit gt 0 or dropped) then
        // An unclassified (still-unsplit product) or CurveFromInvariants-rejected quotient:
        // the analytic side did not fully account for its survivors, so certification is
        // unavailable regardless of how the counts compare.
        vprintf Gluing: "gluingCompositeCRT: unclassified or rejected quotient(s) at level %o, certificate withheld\n", n;
        error if strict,
            Sprintf("Proof := true cannot certify n=%o: the level-n recognition left a quotient unclassified or reconstruction-rejected; use Proof := \"Auto\" or raise Precision", n);
        allCertified := false;
    end if;

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
    strictGuard(strict, info);
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
completeness certificate (exact.m, GaloisStableGluings). "Auto" attempts the supported
exact-count comparisons (prime blocks up to ell = 13) and otherwise returns "traces-only";
false skips them. true is STRICT: it requires non-traces-only evidence for every requested
block (congruence-proved emptiness, reduction to a proved-empty prime block, or an exact
stable-quotient count matched by the analytic count) and raises rather than returning a
traces-only result, hence it errors when the exact layer declines, when the algorithm skips
the comparison (n = 2 under "Periods"), on a Q- or geometrically isogenous pair, or for a
productive e >= 2 or composite block with no v1 proof. The certificate's scope is the paper's
Algorithm 3.4 input contract: geometrically NONisogenous pairs. A pair with IsIsogenous or
GeometricallyIsogenous true (the "isogenous to a square" case) is never certified. Under
"Auto" its proof stays "traces-only" and, at prime blocks where the certificate would have
run, the exact Galois-stable count is still recorded in the block tuple as diagnostic data
(certified false, never compared); under strict true it raises instead. The exact layer runs
at each block's PRIME level only; for e >= 2 blocks it is out of scope in v1 (division
polynomials at ell^e are too large), so those blocks are "traces-only" (raising under strict
true) EXCEPT that a prime level certified empty forces ell^e empty by reduction (certified). The same shortcut applies to a composite n (gluingCompositeCRT): one
certified-empty block (congruence-obstructed or exact-certified-empty) proves the whole
composite empty, and the composite returns "certified" immediately with the remaining
blocks in the tuple reported as unexamined sentinels, not independently certified. The
GluingInfo fields are n, proof ("certified" when every block certified,
else "traces-only"), blocks (one <ell, e, stable_count, analytic_count, certified> tuple
per prime-power block; stable_count is the exact Galois-stable QUOTIENT count when
certified, the DIAGNOSTIC exact count when the pair is geometrically isogenous and the
certificate would otherwise have run (certified false; -1 if the exact layer declined),
and -1 otherwise (including every non-empty e >= 2 block); analytic_count the
block's rational-looking survivor orbit count; the composite certificate compares the
PRODUCT of block graph counts against the analytic graph count and a disagreement is a hard
error), psis (a gluing matrix per returned curve; empty on the Algebraic/BHLS dispatch, which
supplies no matrices), products (recognized <j1, j2> pairs of
the product-type quotients; reported only, never part of any certificate count),
precision (digits of the successful analytic pass), and
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
    // Scope (paper Algorithm 3.4): the certificate is only run for GEOMETRICALLY
    // nonisogenous E1, E2, i.e. Hom_Qbar(E1, E2) = 0. When Hom_Qbar(E1, E2) != 0
    // (a Q-isogeny, or the twist-family/CM degeneracies GeometricallyIsogenous
    // detects, e.g. 54a1 x 54b1 with equal j = 9261/8) the surface E1 x E2
    // carries extra endomorphisms: Aut(E1 x E2) can exceed {+-1}^2, Kani's
    // clean-locus vacuity for product quotients fails (a stable graph CAN have a
    // product-type quotient, e.g. the twisted Weil restriction at n = 2), and
    // distinct M/-M orbits can collide in moduli. The analytic rational-quotient
    // count is then not comparable to the exact stable count, so excluded pairs
    // stay "traces-only"; where the certificate WOULD have run, the exact stable
    // count is still recorded in the block tuple as diagnostics (certified false,
    // never compared, never an error). This is the paper's still-in-development
    // "isogenous to a square" case.
    // Over a number field the exact completeness layer (exact.m) and BHLS both
    // require Q, and IsIsogenous/GeometricallyIsogenous are Q-only, so force the
    // certificate off: number-field blocks are traces-only in v1 (the congruence
    // certified-empty path below still applies, GluingModulus holding over K).
    // gateExact short-circuits, so neither predicate is evaluated over K or when
    // the Proof gates already fail.
    gateExact := (not isNF)
        and ((Proof cmpeq true) or (Proof cmpeq "Auto" and ell le 13));
    excluded := gateExact and (IsIsogenous(E1, E2) or GeometricallyIsogenous(E1, E2));
    runExact := gateExact and not excluded;
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

        // Strict Proof := true (option (a)): only the certified-empty-by-reduction above can
        // satisfy strict at e >= 2. A productive prime power n = ell^e with e >= 2 has no exact
        // layer over the prime level in v1, so everything past this point finishes traces-only
        // (the empty-lift-sweep return and the final e >= 2 return below). Raise here, before
        // the sweep, rather than compute a result strict will not accept.
        error if strict and excluded,
            Sprintf("Proof := true cannot certify n=%o=%o^%o: E1 and E2 are Q-isogenous or geometrically isogenous (outside the certificate's scope); use Proof := \"Auto\"", n, ell, e);
        error if strict,
            Sprintf("Proof := true cannot certify a prime power n=%o=%o^%o with e >= 2: v1 has no exact layer above the prime level (only a prime-level certified-empty reduces to certified); use Proof := \"Auto\"", n, ell, e);

        // The sweep evaluates quotients up to the top level n = ell^e, so its precision is
        // keyed to the top level (honoring a user Precision), not the base prime ell;
        // recognition then reruns on the survivors at that same target precision.
        sweepPrec := (Precision cmpeq false) select GluingPrecisionHeuristic(E1, E2, n) else Precision;
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
            strictGuard(strict, info);   // unreachable under strict (early guard above), belt-and-braces
            return emptyCurves(), info;
        end if;
        prec := (Precision cmpeq false) select GluingPrecisionHeuristic(E1, E2, n) else Precision;
        recIC, recPsi, products, usedPrec := gluingRecognizeAtLevel(E1, E2, survivors, n, prec);
        cs, psisOut, analyticCount := gluingEmitCurves(recIC, recPsi, E1, E2, TraceBound);
        // No exact ell^e layer in v1, so stable_count is -1 and the block is traces-only
        // even when curves are emitted (each still passes the Euler-factor twist
        // certificate to TraceBound in CurveFromInvariants). Recognized products are
        // reported in info`products but never counted (file header).
        info := rec< gluingInfoFmt() | n := n, proof := "traces-only",
            blocks := [<ell, e, -1, analyticCount, false>], psis := psisOut, products := products,
            precision := usedPrec, tracebound := TraceBound >;
        strictGuard(strict, info);   // unreachable under strict (early guard above), belt-and-braces
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
        // BHLS is a complete algorithm for the (n, n)-jacobian gluings at n in {2, 3},
        // and on the certificate's clean locus every stable quotient is a jacobian
        // (product-type graph quotients are impossible there by Kani, file header), so
        // #cs is the analytic count the exact quotient count must match. Jacobian
        // recognition from the period enumeration is deliberately NOT used for the
        // count here (at small mod-n image the conjugation filter barely cuts the pool
        // and look-alike quotients isogenous to E1 x E2 twist-validate); the
        // enumeration runs only to report the product-type quotients' j-pairs in
        // info`products (user-facing data on clean pairs, square-case data on excluded
        // ones), under the same gates the certificate/diagnostics obey.
        products := []; usedPrec := 0; unsplit := 0;
        if runExact or excluded then
            prec := (Precision cmpeq false) select GluingPrecisionHeuristic(E1, E2, n) else Precision;
            _, _, products, usedPrec, unsplit := gluingAnalyticEnumerate(E1, E2, n, prec);
        end if;
        analyticCount := #cs;
        if runExact and unsplit eq 0 then
            block, blockProof := GluingCertificateBlock(E1, E2, n, analyticCount, strict);
        else
            // Cannot certify. Under strict Proof := true (option (a)) that is a hard error
            // with a site-specific message; otherwise record the diagnostic exact count
            // (excluded pair, never compared) or withhold on a residual unsplit product, and
            // finish traces-only. (n = 2 is fine here: the BHLS list #cs IS the certifiable
            // count, so this branch is reached only when excluded or a product is unsplit.)
            error if strict and excluded,
                Sprintf("Proof := true cannot certify n=%o: E1 and E2 are Q-isogenous or geometrically isogenous (the paper's isogenous-to-a-square case, outside the certificate's scope); use Proof := \"Auto\"", n);
            error if strict,
                Sprintf("Proof := true cannot certify n=%o: the analytic enumeration left a product-type quotient unclassified; use Proof := \"Auto\" or raise Precision", n);
            stable := -1;
            if excluded then stable := GaloisStableGluings(E1, E2, n); end if;
            block := <n, 1, stable, analyticCount, false>; blockProof := "traces-only";
        end if;
        info := rec< gluingInfoFmt() | n := n, proof := blockProof,
            blocks := [block], psis := [], products := products,
            precision := usedPrec, tracebound := TraceBound >;
        strictGuard(strict, info);
        return cs, info;
    end if;

    // Analytic (periods) path. The two-pass precision sweep enumerating and
    // recognizing the anti-symplectic quotients is gluingAnalyticEnumerate above;
    // here we reconstruct Q-models from the recognized jacobian invariants and run
    // the completeness certificate.
    prec := (Precision cmpeq false) select GluingPrecisionHeuristic(E1, E2, n) else Precision;
    recIC, recPsi, products, usedPrec, unsplit := gluingAnalyticEnumerate(E1, E2, n, prec);

    // Phase 2 (reconstruction + count), shared with the prime-power path via
    // gluingEmitCurves: one pass over the distinct recognized jacobian moduli. A
    // twist-validated quotient is emitted as a curve and counted; a quotient
    // CurveFromInvariants rejects is not counted and makes the block uncertifiable
    // (dropped flag). Products are reported in info`products but never counted
    // (clean-locus Kani vacuity, file header). analyticCount is the distinct emitted
    // jacobian moduli, the same unit as the exact layer's M/-M orbits.
    cs, psisOut, analyticCount, dropped := gluingEmitCurves(recIC, recPsi, E1, E2, TraceBound);

    // The analytic recognition is a reliable rational-quotient count only at odd n,
    // where the conjugation filter is tight; at n = 2 it barely cuts the candidate
    // pool and look-alike jacobians isogenous to E1 x E2 twist-validate, so the n = 2
    // certificate is left to the BHLS path (which Auto always takes at n = 2, non-
    // isomorphic). Under Algorithm := "Periods" an n = 2 block stays "traces-only".
    // The certificate also requires every survivor accounted for: a residual unsplit
    // product or a CurveFromInvariants-rejected quotient withholds it (traces-only).
    if runExact and n ge 3 and unsplit eq 0 and not dropped then
        block, blockProof := GluingCertificateBlock(E1, E2, n, analyticCount, strict);
    else
        // Cannot certify this block on the analytic path. Under strict Proof := true
        // (option (a)) that is a hard error with a site-specific message; otherwise record
        // the exact stable count as diagnostics (excluded pair, never compared) and finish
        // traces-only. The three strict cases: an excluded (Q-/geometrically isogenous) pair;
        // n = 2, where the analytic recognition is not a reliable count and the certificate
        // is deliberately left to the BHLS path (so "Periods" declines it); or a residual
        // unsplit/reconstruction-rejected quotient at n >= 3.
        error if strict and excluded,
            Sprintf("Proof := true cannot certify n=%o: E1 and E2 are Q-isogenous or geometrically isogenous (the paper's isogenous-to-a-square case, outside the certificate's scope); use Proof := \"Auto\"", n);
        error if strict and n eq 2,
            "Proof := true cannot certify n = 2 on the analytic (Periods) path, which deliberately skips the exact comparison; use Algorithm := \"Auto\" or \"Algebraic\", or Proof := \"Auto\"";
        error if strict,
            Sprintf("Proof := true cannot certify n=%o: the analytic recognition left a quotient unclassified or reconstruction-rejected; use Proof := \"Auto\" or raise Precision", n);
        stable := -1;
        if excluded and n ge 3 then stable := GaloisStableGluings(E1, E2, n); end if;
        block := <n, 1, stable, analyticCount, false>; blockProof := "traces-only";
    end if;
    info := rec< gluingInfoFmt() | n := n, proof := blockProof,
        blocks := [block], psis := psisOut, products := products,
        precision := usedPrec, tracebound := TraceBound >;
    strictGuard(strict, info);
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
