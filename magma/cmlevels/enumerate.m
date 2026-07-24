/*
 * Candidate-ideal enumeration and local-principality filter for the CM
 * gluing-levels engine (magma/cmlevels).
 *
 * A candidate is an integral left O-ideal a with reduced norm n:
 *     n*O <= a <= O,   [O:a] = n^2.
 * It is carried around as its 4x4 integral coordinate matrix X w.r.t. Basis(O)
 * (rows = a Z-basis of a in O-coordinates, Hermite-normalized), together with O.
 *
 * This file is authored across four tasks that all EXTEND it:
 *   T2c (this section) : the local-principality layer -- premise checks, the
 *                        exact gcd test nu(a), and the per-prime cyclicity
 *                        witnesses, with the gcd-vs-cyclicity consistency assert.
 *   T2a  : good-prime catalogues (p+1 maximal-subspace descent).
 *   T2b  : bad-prime bounded descent.
 *   T2d  : CRT assembly + completeness asserts.
 * Every enumerated candidate from T2a/T2b is filtered through the T2c layer, so
 * the helpers below are exported with docstrings for those tasks to call.
 *
 * Mathematics (consensus of gpt-backend-1 sec 2-3, gpt-backend-2 sec 2,
 * gpt-magma-impl sec 3; PLAN section 3):
 *  - Local principality GCD test. With Z-basis b_1..b_4 of a,
 *      nu(a) = gcd( {nrd(b_i)} union {nrd(b_i+b_j) - nrd(b_i) - nrd(b_j) : i<j} ).
 *    For p^e || n, right multiplication by any x in a has determinant nrd(x)^2 and
 *    O_p x <= a_p, forcing v_p(nrd x) >= e; equality is attainable iff a_p is
 *    principal, i.e. iff v_p(nu) = e. Hence a is locally principal everywhere iff
 *    nu(a) = n. No Gorenstein hypothesis. The cross-term form (never the trace-Gram,
 *    whose diagonal is 2*nrd) keeps the gcd exact at p = 2.
 *  - Cyclicity witness. a_p is principal over O_p iff a/pa is cyclic over O/pO
 *    (Nakayama; p in rad(O_p), a full-rank generator is invertible in B_p). With
 *    A_i the matrix of left multiplication by the i-th O-basis element written in
 *    the basis b_1..b_4 of a (integral exactly by left stability), u = sum v_i b_i
 *    generates a/pa iff det[ v*A_1 | ... | v*A_4 ] != 0 in F_p. That determinant is
 *    homogeneous of degree 4 with degree <= 4 in each variable, so vanishing on the
 *    grid {0,1,2,3,4}^4 proves it is the zero polynomial: the grid search is a
 *    COMPLETE decision procedure for p >= 5, and all p^4 vectors are exhausted for
 *    p in {2,3}.
 */

// ---------------------------------------------------------------------------
// Filter-result record schema (T2a/T2b attach these to catalogue entries).
//
// CMPrimeVerdictFormat  (one per prime p | n):
//   p                : the prime.
//   e                : v_p(n) = the exponent p^e || n.
//   LocallyPrincipal : boolean, a_p principal over O_p (the gcd verdict; when
//                      witnesses are requested it equals the cyclicity verdict,
//                      enforced by a hard consistency assert).
//   Searched         : boolean, whether the cyclicity search ran for this prime
//                      (true iff LocalPrincipalityData had Witnesses := true).
//                      Witness/WitnessLift are meaningful ONLY when Searched; when
//                      false they hold [] and 0 as placeholders, which is NOT a
//                      "no witness / nonprincipal" verdict (read LocallyPrincipal).
//   Witness          : the F_p generator vector v as a SeqEnum of integers in
//                      [0,p-1] when principal, or [] ("none") when not; meaningful
//                      only when Searched. A LOCAL generation claim: u = sum v_i b_i
//                      generates a_p, never a global generator.
//   WitnessLift      : the algebra element u = sum v_i b_i when principal, else 0;
//                      meaningful only when Searched.
//
// CMFilterFormat  (one per candidate a):
//   n                : reduced norm, sqrt([O:a]).
//   nu               : nu(a) from the ten integers; a multiple of n, equal iff
//                      locally principal everywhere.
//   LocallyPrincipal : boolean, nu eq n.
//   PerPrime         : List of CMPrimeVerdictFormat records, one per p | n.
// ---------------------------------------------------------------------------
CMPrimeVerdictFormat := recformat< p, e, LocallyPrincipal, Searched, Witness, WitnessLift >;
CMFilterFormat       := recformat< n, nu, LocallyPrincipal, PerPrime >;

// ===========================================================================
// LOCAL-PRINCIPALITY LAYER (T2c)
// ===========================================================================

// ---------------------------------------------------------------------------
// Internal helpers (file-local).
// ---------------------------------------------------------------------------

// Multiplication data of O, as pure functions of O. Returns:
//   B     : the algebra Algebra(O),
//   E     : Basis(O) coerced into B (a Z-basis of O as algebra elements),
//   BO    : rows = coordinates of E in B's standard basis,
//   BOinv : BO^-1 (coordinatizes any element of B against E),
//   AO    : the four row-convention left-multiplication matrices of E on O
//           (AO[i] has row l = O-coordinates of E[i]*E[l]); integral because O is
//           an order, which this asserts.
function OrderMultData(O)
    B := Algebra(O);
    E := [B ! e : e in Basis(O)];
    BO := Matrix(Rationals(), [Eltseq(e) : e in E]);
    BOinv := BO^-1;
    AO := [];
    for i in [1..4] do
        P := Matrix(Rationals(), [Eltseq(E[i]*E[l]) : l in [1..4]]);
        Ai := P * BOinv;
        assert &and[Denominator(c) eq 1 : c in Eltseq(Ai)];
        Append(~AO, Matrix(Integers(), 4, 4, [Integers()!c : c in Eltseq(Ai)]));
    end for;
    return B, E, BO, BOinv, AO;
end function;

// Z-basis of a as algebra elements: row k of Xi is sum_l Xi[k,l]*E[l].
function IdealBasisElts(E, Xi)
    return [ &+[Xi[k,l]*E[l] : l in [1..4]] : k in [1..4] ];
end function;

// nu(a) = gcd of the four reduced norms and the six cross terms
// nrd(b_i+b_j) - nrd(b_i) - nrd(b_j).  Computed exactly (no trace-Gram, no /2).
function NuFromBasis(bk)
    ns := [ Integers() ! Norm(bk[k]) : k in [1..4] ];
    cs := [ Integers() ! (Norm(bk[a]+bk[b]) - Norm(bk[a]) - Norm(bk[b]))
            : a, b in [1..4] | a lt b ];
    return GCD(ns cat cs);
end function;

// Minimal input validation shared by every entry point: X is a 4x4 integral
// full-rank matrix. Returns ok, message, Xi (integer 4x4; zero on failure).
function ValidateIntegral4x4(X)
    Z4 := ZeroMatrix(Integers(), 4, 4);
    if Nrows(X) ne 4 or Ncols(X) ne 4 then
        return false, "X must be a 4x4 matrix", Z4;
    end if;
    if not &and[Denominator(Rationals()!c) eq 1 : c in Eltseq(X)] then
        return false, "X must be integral (containment a subset O: each row must lie in O)", Z4;
    end if;
    Xi := Matrix(Integers(), 4, 4, [Integers()!(Rationals()!c) : c in Eltseq(X)]);
    if Determinant(Xi) eq 0 then
        return false, "X must have full rank 4", Z4;
    end if;
    return true, "", Xi;
end function;

// Left-O-stability: the four matrices Xi*AO[i]*Xi^-1 (left multiplication by the
// i-th O-basis element, in the a-basis) must be integral. Returns ok, message,
// and the integer matrices AA (empty on failure). Their integrality IS left
// stability, and they are the mod-p A_i of the cyclicity test.
function ValidateStable(Xi, AO)
    XiQ := ChangeRing(Xi, Rationals());
    XiInv := XiQ^-1;
    AA := [];
    for i in [1..4] do
        M := XiQ * ChangeRing(AO[i], Rationals()) * XiInv;
        if not &and[Denominator(c) eq 1 : c in Eltseq(M)] then
            return false,
                "a must be left-O-stable (X*mat(e)*X^-1 integral for each O-basis element e)",
                [];
        end if;
        Append(~AA, Matrix(Integers(), 4, 4, [Integers()!c : c in Eltseq(M)]));
    end for;
    return true, "", AA;
end function;

// ---------------------------------------------------------------------------
// Exported: coordinate converters (rows <-> algebra elements).
// ---------------------------------------------------------------------------

intrinsic CMIdealElements(O::AlgQuatOrd, X::Mtrx) -> SeqEnum
{ The Z-basis of the candidate ideal a as algebra elements: row k of the integral
  4x4 matrix X gives sum_l X[k,l]*Basis(O)[l]. Requires X integral, 4x4, full rank. }
    ok, msg, Xi := ValidateIntegral4x4(X);
    require ok: msg;
    _, E := OrderMultData(O);
    return IdealBasisElts(E, Xi);
end intrinsic;

intrinsic CMCoordinateMatrix(O::AlgQuatOrd, elts::SeqEnum) -> Mtrx
{ Integer matrix whose rows are the O-coordinates (w.r.t. Basis(O)) of the given
  elements, which must lie in O. Inverts CMIdealElements on a Z-basis of a. }
    B, _, _, BOinv := OrderMultData(O);
    rows := [];
    for x in elts do
        c := Vector(Rationals(), Eltseq(B!x)) * BOinv;
        require &and[Denominator(t) eq 1 : t in Eltseq(c)]: "each element must lie in O";
        Append(~rows, [Integers()!t : t in Eltseq(c)]);
    end for;
    return Matrix(Integers(), rows);
end intrinsic;

intrinsic CMLeftIdealMatrix(O::AlgQuatOrd, gens::SeqEnum) -> Mtrx
{ The Hermite-normalized 4x4 O-coordinate matrix of the left ideal O*gens generated
  by the given algebra elements (which must lie in O): the Z-generators
  Basis(O)[i]*g are coordinatized, stacked, and Hermite-reduced. Requires the
  result to have full rank 4 (finite index in O). A convenience for building
  candidates; T2a/T2b build X directly from submodule enumeration. }
    B, E, _, BOinv := OrderMultData(O);
    rows := [];
    for g in gens do
        gg := B ! g;
        for i in [1..4] do
            c := Vector(Rationals(), Eltseq(E[i]*gg)) * BOinv;
            require &and[Denominator(t) eq 1 : t in Eltseq(c)]:
                "O*gens must be integral: each generator must lie in O";
            Append(~rows, [Integers()!t : t in Eltseq(c)]);
        end for;
    end for;
    H := HermiteForm(Matrix(Integers(), rows));
    nz := [ H[l] : l in [1..Nrows(H)] | not IsZero(H[l]) ];
    require #nz eq 4: "O*gens must have full rank 4 (finite index in O)";
    return HermiteForm(Matrix(nz));
end intrinsic;

// ---------------------------------------------------------------------------
// Exported: premise checks (called by every decision below).
// ---------------------------------------------------------------------------

intrinsic CMCandidatePremise(O::AlgQuatOrd, X::Mtrx, n::RngIntElt) -> Mtrx, SeqEnum, SeqEnum
{ Premise checks for a candidate ideal a given by its integral 4x4 coordinate
  matrix X w.r.t. Basis(O) and claimed reduced norm n. Requires, with clean errors:
  X integral (a subset O), 4x4, full rank, in Hermite normal form; [O:a] = det(X)
  = n^2; n*O subset a; and left-O-stability. Returns the coerced integer matrix Xi,
  the Z-basis of a as algebra elements, and the integer left-multiplication matrices
  A_i on a (row l of A_i = a-coordinates of Basis(O)[i]*b_l). Every downstream test
  goes through this first. }
    ok, msg, Xi := ValidateIntegral4x4(X);
    require ok: msg;
    require Xi eq HermiteForm(Xi): "X must be in Hermite normal form (canonical row basis)";
    require n ge 1: "n must be a positive integer";
    require Determinant(Xi) eq n^2: "index [O:a] = det(X) must equal n^2";
    XiQ := ChangeRing(Xi, Rationals());
    require &and[Denominator(c) eq 1 : c in Eltseq(n*XiQ^-1)]:
        "n*O must be contained in a (n*X^-1 must be integral)";
    B, E, _, _, AO := OrderMultData(O);
    ok2, msg2, AA := ValidateStable(Xi, AO);
    require ok2: msg2;
    return Xi, IdealBasisElts(E, Xi), AA;
end intrinsic;

// ---------------------------------------------------------------------------
// Exported: the gcd test and the local-principality decision.
// ---------------------------------------------------------------------------

intrinsic LocalPrincipalityGCD(O::AlgQuatOrd, X::Mtrx) -> RngIntElt
{ nu(a) = gcd of the four reduced norms nrd(b_i) and the six cross terms
  nrd(b_i+b_j)-nrd(b_i)-nrd(b_j), for the Z-basis b_1..b_4 of the candidate a given
  by the integral 4x4 matrix X. A multiple of n = sqrt([O:a]) for a genuine left
  ideal, equal to n iff a is locally principal everywhere. Requires only X integral,
  4x4, full rank (nu is defined for any such lattice). }
    ok, msg, Xi := ValidateIntegral4x4(X);
    require ok: msg;
    _, E := OrderMultData(O);
    nu := NuFromBasis(IdealBasisElts(E, Xi));
    assert nu gt 0;
    return nu;
end intrinsic;

intrinsic IsLocallyPrincipal(O::AlgQuatOrd, X::Mtrx, n::RngIntElt) -> BoolElt
{ True iff the candidate a (integral 4x4 X, claimed reduced norm n) is locally
  principal at every prime, i.e. nu(a) = n. Runs the full premise checks first. }
    _, bk := CMCandidatePremise(O, X, n);
    nu := NuFromBasis(bk);
    assert nu gt 0;
    assert IsDivisibleBy(nu, n);   // norm-index identity forces n | nu post-premise
    return nu eq n;
end intrinsic;

// ---------------------------------------------------------------------------
// Exported: per-prime cyclicity witness.
// ---------------------------------------------------------------------------

intrinsic CyclicityWitness(O::AlgQuatOrd, X::Mtrx, p::RngIntElt) -> BoolElt, ModTupFldElt, AlgQuatElt
{ Decide whether a/pa is cyclic over O/pO (equivalently a_p principal over O_p) for
  the candidate a (integral 4x4 X) at the prime p, by a COMPLETE deterministic
  search for v in F_p^4 with det[ v*A_1 | ... | v*A_4 ] != 0, where A_i is left
  multiplication by Basis(O)[i] on a reduced mod p. Search order: the four standard
  basis vectors, then all-ones, then the full grid using all of 0..p-1 per
  coordinate for p in 2 or 3, or the five values 0,1,2,3,4 per coordinate for
  p >= 5 (a complete test set for the degree-4 determinant). On success returns
  true, the witness vector v over F_p, and its lift u = sum v_i b_i; otherwise
  false, the zero vector, and 0. v is a PER-PRIME local generation claim. Asserts
  the A_i are integral, which re-proves left stability. }
    require IsPrime(p): "p must be prime";
    ok, msg, Xi := ValidateIntegral4x4(X);
    require ok: msg;
    B, E, _, _, AO := OrderMultData(O);
    oks, msgs, AA := ValidateStable(Xi, AO);
    require oks: msgs;
    // Frozen consistency requirement on EVERY path (direct callers included): the
    // cyclicity verdict must agree with the gcd verdict. n comes from the premise
    // index [O:a] = det(X) = n^2 (asserted here since we skip CMCandidatePremise).
    dX := Determinant(Xi);
    n := Isqrt(dX);
    assert n^2 eq dX;
    e := Valuation(n, p);
    nu := LocalPrincipalityGCD(O, X);
    expected := (Valuation(nu, p) eq e);
    F := GF(p);
    AAp := [ ChangeRing(AA[i], F) : i in [1..4] ];
    bk := IdealBasisElts(E, Xi);
    Sset := (p le 3) select [0..p-1] else [0, 1, 2, 3, 4];
    cands := [ [1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1], [1,1,1,1] ]
             cat [ [t[1], t[2], t[3], t[4]] : t in CartesianPower(Sset, 4) ];
    found := false;
    vwit := VectorSpace(F, 4) ! 0;
    uwit := B ! 0;
    for vv in cands do
        v := Vector(F, vv);
        if IsZero(v) then continue; end if;
        D := Matrix(F, [ Eltseq(v*AAp[i]) : i in [1..4] ]);
        if not IsZero(Determinant(D)) then
            found := true;
            vwit := v;
            uwit := &+[ (Integers()!v[l]) * bk[l] : l in [1..4] ];
            break;
        end if;
    end for;
    error if found ne expected,
        "CONSISTENCY FAILURE at p=" cat IntegerToString(p) cat
        ": gcd test (nu=" cat IntegerToString(nu) cat ", n=" cat IntegerToString(n) cat
        ") says " cat (expected select "principal" else "not principal") cat
        " but cyclicity search found " cat (found select "a witness" else "none");
    return found, vwit, uwit;
end intrinsic;

// ---------------------------------------------------------------------------
// Exported: the full filter record, with the consistency cross-check.
// ---------------------------------------------------------------------------

intrinsic LocalPrincipalityData(O::AlgQuatOrd, X::Mtrx, n::RngIntElt : Witnesses := true) -> Rec
{ The filter result for the candidate a (integral 4x4 X, reduced norm n): a
  CMFilterFormat record with n, nu(a), the global verdict nu eq n, and a per-prime
  List of CMPrimeVerdictFormat records over p | n. When Witnesses is true (default),
  a cyclicity search runs at each p | n and its verdict is cross-checked against the
  gcd verdict; a disagreement is a hard error (mathematical inconsistency, never
  swallowed). Runs the full premise checks first. }
    Xi, bk, AA := CMCandidatePremise(O, X, n);
    B := Algebra(O);
    nu := NuFromBasis(bk);
    assert nu gt 0;
    assert IsDivisibleBy(nu, n);   // norm-index identity forces n | nu post-premise
    lp := (nu eq n);

    perprime := [* *];
    allp := true;
    for p in PrimeDivisors(n) do
        ep := Valuation(n, p);
        gcdverdict := (Valuation(nu, p) eq ep);
        allp := allp and gcdverdict;
        wit := [ Integers() | ];
        witlift := B ! 0;
        if Witnesses then
            found, v, u := CyclicityWitness(O, X, p);
            error if found ne gcdverdict,
                "CONSISTENCY FAILURE at p=" cat IntegerToString(p) cat
                ": gcd test says " cat (gcdverdict select "principal" else "not principal") cat
                " but cyclicity search found " cat (found select "a witness" else "none");
            if found then
                // A generator of a/pa lifts to a local generator O_p u = a_p, so
                // its reduced norm has p-valuation exactly e_p; cross-check that.
                assert Valuation(Integers()!Norm(u), p) eq ep;
                wit := [ Integers()!v[l] : l in [1..4] ];
                witlift := u;
            end if;
        end if;
        Append(~perprime, rec< CMPrimeVerdictFormat |
            p := p, e := ep, LocallyPrincipal := gcdverdict, Searched := Witnesses,
            Witness := wit, WitnessLift := witlift >);
    end for;
    assert lp eq allp;   // nu eq n iff locally principal at every p | n

    return rec< CMFilterFormat |
        n := n, nu := nu, LocallyPrincipal := lp, PerPrime := perprime >;
end intrinsic;
