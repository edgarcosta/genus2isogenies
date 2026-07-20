/*
 * From numeric invariants back to a curve over Q.
 *
 * The analytic engine (analytic.m) hands us, for each anti-symplectic gluing,
 * a numeric Igusa-Clebsch quadruple. Two steps recover an exact curve:
 *
 *   1. RecognizeIgusaClebsch: read the quadruple as a point of the weighted
 *      projective space P(2, 4, 6, 10) and recognize a rational representative.
 *      We pivot on I2 (weight 2 divides 4, 6 and 10, so the normalization uses
 *      only even powers of a single well-defined scalar 1/I2 and carries no
 *      sign/root ambiguity). The other coordinates are the unique clean pivot;
 *      a quotient with I2 projectively zero (|I2^5/I10| below the gate) is not
 *      handled here and is reported unrecognized.
 *
 *   2. CurveFromInvariants: reconstruct one Q-model with
 *      HyperellipticCurveFromIgusaClebsch (correct up to a quadratic twist),
 *      then pin the twist by matching Euler factors against E1 x E2. The
 *      surviving twist is certificate (a): its L-factor equals L(E1) L(E2) at
 *      every good prime up to TraceBound.
 */

// I2 pivot in P(2,4,6,10). Returns false when I10 vanishes (not a Jacobian) or
// I2 is projectively zero (|I2^5/I10| below the recognition gate), else the
// numeric normalized quadruple [1, I4/I2^2, I6/I2^3, I10/I2^5].
function normalizeIC(IC)
    I2 := IC[1]; I4 := IC[2]; I6 := IC[3]; I10 := IC[4];
    CC := Parent(I2);
    P := Floor(Precision(CC));
    gate := 10^(-(P div 2));
    if Abs(I10) le gate or Abs(I2^5 / I10) lt gate then
        return false, [CC | ];
    end if;
    return true, [CC | 1, I4 / I2^2, I6 / I2^3, I10 / I2^5];
end function;

intrinsic RecognizeIgusaClebsch(IC::SeqEnum) -> BoolElt, SeqEnum
{Recognize the numeric Igusa-Clebsch quadruple IC = [I2, I4, I6, I10] as a
rational point of the weighted projective space P(2, 4, 6, 10). Normalizes by
the I2 pivot (the only weight dividing all of 4, 6, 10) and RecognizeRational's
each coordinate, all-or-nothing. On success returns true and a rational quadruple
[1, I4/I2^2, I6/I2^3, I10/I2^5] with nonzero last coordinate (a bona fide
Igusa-Clebsch quadruple for HyperellipticCurveFromIgusaClebsch); otherwise false.}
    require #IC eq 4: "IC must be a length-4 Igusa-Clebsch quadruple";
    ok, Q := normalizeIC(IC);
    if not ok then return false, [Rationals() | ]; end if;
    out := [Rationals() | ];
    for q in Q do
        okq, r := RecognizeRational(q);
        if not okq then return false, [Rationals() | ]; end if;
        Append(~out, r);
    end for;
    if out[4] eq 0 then return false, [Rationals() | ]; end if;
    return true, out;
end intrinsic;

intrinsic IgusaClebschNearRational(IC::SeqEnum) -> BoolElt
{True when the I2-normalized coordinates of the numeric quadruple IC all have
imaginary part below the recognition gate 10^(-P/2) (P the working precision):
the quotient looks defined over Q, so a RecognizeIgusaClebsch failure is a
precision shortfall worth a retry rather than a genuinely complex (conjugate)
quotient. False when I2 is projectively zero or a coordinate is visibly complex.}
    require #IC eq 4: "IC must be a length-4 Igusa-Clebsch quadruple";
    ok, Q := normalizeIC(IC);
    if not ok then return false; end if;
    P := Floor(Precision(Universe(Q)));
    gate := 10^(-(P div 2));
    return forall{q : q in Q | Abs(Im(q)) le gate};
end intrinsic;

intrinsic CurveFromInvariants(ICQ::SeqEnum, E1::CrvEll, E2::CrvEll : TraceBound := 1000) -> BoolElt, CrvHyp
{Reconstruct a genus-2 curve over Q from the rational Igusa-Clebsch quadruple ICQ
and pin its quadratic twist against E1 x E2. C0 := HyperellipticCurveFromIgusaClebsch(ICQ)
is correct up to a quadratic twist; the twisting discriminant is searched over
squarefree d (both signs) supported on the primes dividing 2, cond(E1), cond(E2).
A d survives when EulerFactor(QuadraticTwist(C0, d), p) = EulerFactor(E1, p)
EulerFactor(E2, p) at the first 5 good primes; the survivors are then checked at
every good p <= TraceBound (certificate (a)). Returns false when nothing survives
(the recognition was spurious); errors when two non-isomorphic curves survive
(ambiguous twist).}
    require #ICQ eq 4: "ICQ must be a length-4 Igusa-Clebsch quadruple";
    require Type(BaseRing(E1)) eq FldRat and Type(BaseRing(E2)) eq FldRat:
        "E1 and E2 must be defined over Q";
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

    // Candidate discriminants: +-1 times squarefree products of the bad primes
    // of the pair, of the reconstructed model, and 2. The connecting twist is
    // supported there: Jac(C) is isogenous to E1 x E2 (conductor over the E_i
    // bad primes), and C0 may carry spurious bad primes from the reconstruction
    // that the twist must cancel.
    badset := {2} join Seqset(BadPrimes(C0))
                    join Seqset(BadPrimes(E1)) join Seqset(BadPrimes(E2));
    primeset := Sort(Setseq(badset));
    ds := [Integers() | 1];
    for p in primeset do ds cat:= [d * p : d in ds]; end for;
    ds cat:= [-d : d in ds];

    // A prime is good for every candidate twist and both E_i once it avoids
    // badset (every candidate d is supported there, so p does not divide d and
    // KroneckerSymbol(d, p) = +-1).

    // For a good p, EulerFactor(QuadraticTwist(C0, d), p) = EulerFactor(C0, p)
    // with T |-> KroneckerSymbol(d, p) * T, so EulerFactor(C0, p) is computed
    // once per prime and reused across all candidate d.
    screen := [];
    p := 1;
    while #screen lt 5 do
        p := NextPrime(p);
        if p notin badset then Append(~screen, p); end if;
    end while;
    Pol := Parent(EulerFactor(C0, screen[1]));
    T := Pol.1;
    efScreen := [EulerFactor(C0, p) : p in screen];
    tgScreen := [EulerFactor(E1, p) * EulerFactor(E2, p) : p in screen];
    survivors := [d : d in ds |
        forall{i : i in [1 .. #screen] |
            Evaluate(efScreen[i], KroneckerSymbol(d, screen[i]) * T) eq tgScreen[i]}];
    if #survivors eq 0 then
        vprintf Gluing: "CurveFromInvariants: no twist survives the 5-prime screen (spurious recognition)\n";
        return false, C0;
    end if;

    // Full trace certificate at every good p <= TraceBound.
    fullp := [p : p in PrimesUpTo(TraceBound) | p notin badset];
    efFull := [EulerFactor(C0, p) : p in fullp];
    tgFull := [EulerFactor(E1, p) * EulerFactor(E2, p) : p in fullp];
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
    // ambiguity.
    reps := [];
    for d in survfull do
        C := QuadraticTwist(C0, d);
        if forall{D : D in reps | not IsIsomorphic(C, D)} then Append(~reps, C); end if;
    end for;
    require #reps eq 1:
        "ambiguous twist: distinct curves match L(E1) L(E2) up to TraceBound";
    return true, reps[1];
end intrinsic;
