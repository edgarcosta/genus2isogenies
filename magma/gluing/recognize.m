/*
 * From numeric invariants back to a curve over Q.
 *
 * The analytic engine (analytic.m) hands us, for each anti-symplectic gluing,
 * a numeric Igusa-Clebsch quadruple. Two steps recover an exact curve:
 *
 *   1. RecognizeIgusaClebsch: read the quadruple as a point of the weighted
 *      projective space P(2, 4, 6, 10) and recognize a rational representative.
 *      A dimensionless projective scale S = max_w |I_w|^(1/w) turns each
 *      coordinate into a scale-free ratio rho_w = |I_w|/S^w whose largest
 *      surviving weight picks a weight-zero chart: the generic I2 pivot
 *      (I4/I2^2, I6/I2^3, I10/I2^5) when I2 is projectively nonzero, else the
 *      ratio chart t1 = I4 I6/I10, t2 = I4^5/I10^2 on the I2 = 0 locus (two
 *      twist-equivalent signs [0, t2, +-t1 t2, +-t2^2]), its I4 = 0 stratum
 *      [0, 0, u^2, u^3] (u = I6^5/I10^3), and its all-zero point [0, 0, 0, 1].
 *      A descent into a lower stratum is gated on the relevant weight-zero ratio
 *      recognizing as EXACT zero, never a magnitude alone. See
 *      docs/superpowers/ws-b-gpt-verdict.md Part C.
 *
 *   2. CurveFromInvariants: reconstruct one Q-model with
 *      HyperellipticCurveFromIgusaClebsch (correct up to a quadratic twist),
 *      then pin the twist by matching Euler factors against E1 x E2. The
 *      surviving twist is certificate (a): its L-factor equals L(E1) L(E2) at
 *      every good prime up to TraceBound.
 */

// RecognizeRational at an explicit precision P (the working decimal precision).
// Magma complex fields carry one field-level precision, not a per-element one, so
// passing Digits := P pins the recognition to the working precision P directly.
function recQ(z, P)
    return RecognizeRational(z : Digits := P);
end function;

// The dimensionless multichart dispatcher shared by RecognizeIgusaClebsch and
// IgusaClebschNearRational. The projective scale S = max_w |I_w|^(1/w) sends
// S |-> |lambda| S under weighted scaling, so rho_w = |I_w|/S^w is scale-free and
// the gate 10^(-P/2) is projectively meaningful. Returns false for a non-Jacobian
// (I10 projectively zero); otherwise the chart name and the complex weight-zero
// chart coordinates a near-real test inspects and recognition reads:
//   "I2"   -> [I4/I2^2, I6/I2^3, I10/I2^5]   (I2 projectively nonzero: generic)
//   "I4"   -> [t1, t2] = [I4 I6/I10, I4^5/I10^2]        (I2 projectively zero)
//   "I6"   -> [u] = [I6^5/I10^3]                     (I2, I4 projectively zero)
//   "zero" -> [] (the single point [0, 0, 0, 1])  (I2, I4, I6 projectively zero)
function icChart(IC)
    I2 := IC[1]; I4 := IC[2]; I6 := IC[3]; I10 := IC[4];
    CC := Parent(I2);
    P := Floor(Precision(CC));
    gate := 10^(-(P div 2));
    S := Max([Abs(I2)^(1/2), Abs(I4)^(1/4), Abs(I6)^(1/6), Abs(I10)^(1/10)]);
    if S eq 0 then return false, "", [CC | ]; end if;
    if Abs(I10) / S^10 le gate then return false, "", [CC | ]; end if;
    if Abs(I2) / S^2 gt gate then
        return true, "I2", [CC | I4 / I2^2, I6 / I2^3, I10 / I2^5];
    end if;
    if Abs(I4) / S^4 gt gate then
        return true, "I4", [CC | I4 * I6 / I10, I4^5 / I10^2];
    end if;
    if Abs(I6) / S^6 gt gate then
        return true, "I6", [CC | I6^5 / I10^3];
    end if;
    return true, "zero", [CC | ];
end function;

// Recognize IC as rational representative(s) of its weighted-projective point.
// The generic I2 chart yields a single candidate [1, I4/I2^2, I6/I2^3, I10/I2^5];
// the I2 = 0 ratio chart yields the two signs [0, t2, +-t1 t2, +-t2^2] (equivalent
// under the weighted scaling lambda = i). Each descent into a lower stratum first
// requires the relevant weight-zero ratio to recognize as EXACT zero, so I2 = 0 is
// never inferred from a magnitude gate alone (verdict Part C).
function recognizeICCandidates(IC)
    P := Floor(Precision(Parent(IC[1])));
    ok, chart, coords := icChart(IC);
    if not ok then return false, []; end if;

    if chart eq "I2" then
        out := [Rationals() | 1];
        for z in coords do
            okz, q := recQ(z, P);
            if not okz then return false, []; end if;
            Append(~out, q);
        end for;
        if out[4] eq 0 then return false, []; end if;
        return true, [out];
    end if;

    // I2 projectively zero: confirm it is exactly zero before leaving the pivot.
    okz2, z2 := recQ(IC[1]^5 / IC[4], P);
    if not okz2 or z2 ne 0 then return false, []; end if;

    if chart eq "I4" then
        ok1, t1 := recQ(coords[1], P);
        ok2, t2 := recQ(coords[2], P);
        if not ok1 or not ok2 or t2 eq 0 then return false, []; end if;
        return true, [ [Rationals() | 0, t2,  t1 * t2,  t2^2],
                       [Rationals() | 0, t2, -t1 * t2, -t2^2] ];
    end if;

    // I4 projectively zero: confirm exact zero before the I6 stratum.
    okz4, z4 := recQ(IC[2]^5 / IC[4]^2, P);
    if not okz4 or z4 ne 0 then return false, []; end if;

    if chart eq "I6" then
        oku, u := recQ(coords[1], P);
        if not oku or u eq 0 then return false, []; end if;
        return true, [ [Rationals() | 0, 0, u^2, u^3] ];
    end if;

    // I6 projectively zero: confirm exact zero, leaving the single point [0,0,0,1].
    okz6, z6 := recQ(IC[3]^5 / IC[4]^3, P);
    if not okz6 or z6 ne 0 then return false, []; end if;
    return true, [ [Rationals() | 0, 0, 0, 1] ];
end function;

intrinsic RecognizeIgusaClebsch(IC::SeqEnum) -> BoolElt, SeqEnum
{Recognize the numeric Igusa-Clebsch quadruple IC = [I2, I4, I6, I10] as a
rational point of the weighted projective space P(2, 4, 6, 10). A dimensionless
projective scale picks a weight-zero chart (the generic I2 pivot when I2 is
projectively nonzero, else the t1 = I4 I6/I10, t2 = I4^5/I10^2 ratio chart on the
I2 = 0 locus and its I4 = 0 / all-zero sub-strata) and recognizes its coordinates
with RecognizeRational. On success returns true and a rational Igusa-Clebsch quadruple with
nonzero last coordinate; otherwise false. On the I2 = 0 ratio chart the two sign
candidates [0, t2, +-t1 t2, +-t2^2] are weighted-projectively equivalent
(lambda = i), so the first that reconstructs via HyperellipticCurveFromIgusaClebsch
is returned and the twist is pinned downstream by CurveFromInvariants.}
    require #IC eq 4: "IC must be a length-4 Igusa-Clebsch quadruple";
    ok, cands := recognizeICCandidates(IC);
    if not ok then return false, [Rationals() | ]; end if;
    if #cands eq 1 then return true, cands[1]; end if;
    // Two I2 = 0 signs: return the first that builds a curve. Both build in
    // practice (same geometric point, lambda = i apart), so this only tolerates a
    // reconstruction convention; the correct twist is pinned later by CurveFromInvariants.
    for c in cands do
        built := false;
        try
            _ := HyperellipticCurveFromIgusaClebsch(c);
            built := true;
        catch e
            built := false;
        end try;
        if built then return true, c; end if;
    end for;
    return false, [Rationals() | ];
end intrinsic;

intrinsic IgusaClebschNearRational(IC::SeqEnum) -> BoolElt
{True when the numeric quadruple IC looks defined over Q: the weight-zero chart
coordinates chosen by the same multichart dispatcher as RecognizeIgusaClebsch (the
generic I2 pivot when available, else the I2 = 0 ratio chart t1, t2 or the I4 = 0
chart u; the all-zero point is rational) all have imaginary part below the gate
10^(-P/2) (P the working precision). So a RecognizeIgusaClebsch failure here is a
precision shortfall worth a retry rather than a genuinely complex (conjugate)
quotient. False for a non-Jacobian (I10 projectively zero) or a visibly complex
chart coordinate.}
    require #IC eq 4: "IC must be a length-4 Igusa-Clebsch quadruple";
    ok, _, coords := icChart(IC);
    if not ok then return false; end if;
    P := Floor(Precision(Parent(IC[1])));
    gate := 10^(-(P div 2));
    return forall{z : z in coords | Abs(Im(z)) le gate};
end intrinsic;

// Rational prime p usable for the quadratic-twist certificate, and the target
// degree-4 Euler factor EulerFactor(E1,.)*EulerFactor(E2,.) (coerced into Pol) to
// compare against EulerFactor(C0^d, p). Over Q (isNF false) every prime is usable
// and . = p. Over a number field K only primes p with a degree-1 unramified prime
// P | p are usable (the "split" primes, v1 restriction) and . = P: there K_P = Q_p,
// so Frob_P = Frob_p on the l-adic Tate module and L_p(Jac C/Q) = L_P(Jac C_K) =
// L_P(E1) L_P(E2) (base change of the K-isogeny Jac C_K ~ E1 x E2). Any such P above
// p gives the same product, so no P1/P2 assignment is needed. OK is MaximalOrder(K)
// over a number field (unused over Q).
function twistTarget(E1, E2, p, isNF, OK, Pol)
    if not isNF then
        return true, Pol ! (EulerFactor(E1, p) * EulerFactor(E2, p));
    end if;
    for f in Factorization(p * OK) do
        if f[2] eq 1 and Degree(f[1]) eq 1 then
            return true, Pol ! (EulerFactor(E1, f[1]) * EulerFactor(E2, f[1]));
        end if;
    end for;
    return false, Pol ! 1;
end function;

intrinsic CurveFromInvariants(ICQ::SeqEnum, E1::CrvEll, E2::CrvEll : TraceBound := 1000) -> BoolElt, CrvHyp
{Reconstruct a genus-2 curve over Q from the rational Igusa-Clebsch quadruple ICQ
and pin its quadratic twist against E1 x E2. C0 := HyperellipticCurveFromIgusaClebsch(ICQ)
is correct up to a quadratic twist; the twisting discriminant is searched over
squarefree d (both signs) supported on the primes dividing 2, cond(E1), cond(E2).
A d survives when EulerFactor(QuadraticTwist(C0, d), p) = EulerFactor(E1, p)
EulerFactor(E2, p) at the first 5 good primes; the survivors are then checked at
every good p <= TraceBound (certificate (a)). Returns false when nothing survives
(the recognition was spurious); errors when two non-isomorphic curves survive
(ambiguous twist).

E1, E2 over a number field K (experimental, Task 13): the reconstruction C0 and
its twists stay over Q (ICQ is rational), and the twist certificate compares at the
"split" rational primes p only, EulerFactor(E1, P) EulerFactor(E2, P) for a
degree-1 unramified prime P | p (twistTarget). This is a necessary condition family
(a family, since only split p are used) that separates quadratic twists; the loud
uniqueness/ambiguity contract is unchanged.}
    require #ICQ eq 4: "ICQ must be a length-4 Igusa-Clebsch quadruple";
    K := BaseRing(E1);
    require BaseRing(E2) cmpeq K: "E1 and E2 must share a base field";
    require Type(K) eq FldRat or ISA(Type(K), FldNum):
        "E1 and E2 must be defined over Q or a common number field";
    isNF := Type(K) ne FldRat;
    C0 := HyperellipticCurveFromIgusaClebsch(ICQ);
    // HyperellipticCurveFromIgusaClebsch can return a model over a quadratic
    // field when the field of moduli is not a field of definition (Mestre
    // obstruction): a Q-gluing at those invariants is not reachable by this
    // (Igusa-Clebsch, quadratic-twist) reconstruction, so drop the quotient.
    if Type(BaseRing(C0)) ne FldRat then
        vprintf Gluing: "CurveFromInvariants: reconstruction landed off Q (Mestre obstruction), dropping\n";
        return false, C0;
    end if;
    C0 := IntegralModel(C0);

    // Extra geometric automorphisms (order > 2) admit non-quadratic twists: the
    // search below only tries quadratic twists, so a wrong-but-isogenous twist
    // can still pass the trace certificate (observed on bhls2 corpus entry 4,
    // n = 2). Only the Algebraic path (n in {2, 3}) currently reconstructs these
    // quotients correctly, so drop rather than risk certifying the wrong curve.
    if #GeometricAutomorphismGroup(C0) gt 2 then
        vprintf Gluing: "CurveFromInvariants: C0 has extra automorphisms (order %o > 2), quadratic-twist pinning cannot certify the right twist, dropping (the Algebraic path handles n in {2, 3})\n",
            #GeometricAutomorphismGroup(C0);
        return false, C0;
    end if;

    // Candidate discriminants: +-1 times squarefree products of the bad primes
    // of the pair, of the reconstructed model, and 2. The connecting twist is
    // supported there: Jac(C) is isogenous to E1 x E2 (conductor over the E_i
    // bad primes), and C0 may carry spurious bad primes from the reconstruction
    // that the twist must cancel.
    // Rational bad primes: over Q, the bad primes of C0 and both E_i (with 2).
    // Over K, BadPrimes is not defined for CrvEll[FldNum], so use the rational
    // primes below the E_i conductors and the ramified primes of K (disc OK).
    if isNF then
        OK := MaximalOrder(K);
        badE := {2} join Seqset(PrimeDivisors(Integers() ! Norm(Conductor(E1))))
                     join Seqset(PrimeDivisors(Integers() ! Norm(Conductor(E2))))
                     join Seqset(PrimeDivisors(Integers() ! Discriminant(OK)));
    else
        OK := Integers();
        badE := {2} join Seqset(BadPrimes(E1)) join Seqset(BadPrimes(E2));
    end if;
    badset := badE join Seqset(BadPrimes(C0));
    primeset := Sort(Setseq(badset));
    ds := [Integers() | 1];
    for p in primeset do ds cat:= [d * p : d in ds]; end for;
    ds cat:= [-d : d in ds];

    // A usable prime avoids badset (every candidate d is supported there, so p
    // does not divide d and KroneckerSymbol(d, p) = +-1) and, over K, splits to a
    // degree-1 unramified prime P | p (twistTarget); the target degree-4 factor to
    // match EulerFactor(C0^d, p) against is EulerFactor(E1,.) EulerFactor(E2,.),
    // . = p over Q and . = P over K. For a good p, EulerFactor(QuadraticTwist(C0,
    // d), p) = EulerFactor(C0, p) with T |-> KroneckerSymbol(d, p) T, so
    // EulerFactor(C0, p) is computed once per prime and reused across candidate d.
    Pol := PolynomialRing(Integers()); T := Pol.1;
    screen := []; tgScreen := [];
    p := 1;
    while #screen lt 5 do
        p := NextPrime(p);
        if p in badset then continue; end if;
        ok, tg := twistTarget(E1, E2, p, isNF, OK, Pol);
        if ok then Append(~screen, p); Append(~tgScreen, tg); end if;
    end while;
    efScreen := [Pol ! EulerFactor(C0, q) : q in screen];
    survivors := [d : d in ds |
        forall{i : i in [1 .. #screen] |
            Evaluate(efScreen[i], KroneckerSymbol(d, screen[i]) * T) eq tgScreen[i]}];
    if #survivors eq 0 then
        vprintf Gluing: "CurveFromInvariants: no twist survives the 5-prime screen (spurious recognition)\n";
        return false, C0;
    end if;

    // Full trace certificate at every usable good p <= TraceBound.
    fullp := []; tgFull := [];
    for p in PrimesUpTo(TraceBound) do
        if p in badset then continue; end if;
        ok, tg := twistTarget(E1, E2, p, isNF, OK, Pol);
        if ok then Append(~fullp, p); Append(~tgFull, tg); end if;
    end for;
    efFull := [Pol ! EulerFactor(C0, q) : q in fullp];
    survfull := [d : d in survivors |
        forall{i : i in [1 .. #fullp] |
            Evaluate(efFull[i], KroneckerSymbol(d, fullp[i]) * T) eq tgFull[i]}];
    if #survfull eq 0 then
        vprintf Gluing: "CurveFromInvariants: no twist survives the full trace check (spurious recognition)\n";
        return false, C0;
    end if;

    // Distinct twisting discriminants can give isomorphic curves (curves whose
    // Jacobian carries extra automorphisms, e.g. CM, admit several L-equivalent
    // twists that coincide); only genuinely non-isomorphic survivors are an
    // ambiguity. Over a number field K this ambiguity is INTRINSIC and expected on
    // a non-empty gluing: C and its twist by disc(K) are K-isomorphic (disc(K) is a
    // square in K), so both have Jac splitting over K as E1 x E2, and NO K-Euler
    // factor (split or inert) separates them (chi_disc(K)(p) = +1 at split p, and
    // squares to +1 at inert p). The K-gluing data thus pins C over Q only up to
    // disc(K)-twist; v1 keeps the loud contract (surfacing this as an error) rather
    // than silently choosing a descent. The corpus nf pair glues to nothing, so it
    // never reaches here; a non-empty K input is out of v1 scope.
    reps := [];
    for d in survfull do
        C := QuadraticTwist(C0, d);
        if forall{D : D in reps | not IsIsomorphic(C, D)} then Append(~reps, C); end if;
    end for;
    require #reps eq 1:
        Sprintf("ambiguous twist: distinct curves match L(E1) L(E2) up to TraceBound (E1 = %o, E2 = %o)",
            aInvariants(E1), aInvariants(E2));
    return true, reps[1];
end intrinsic;
