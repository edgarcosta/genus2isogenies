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

// ===========================================================================
// GOOD-PRIME CATALOGUES (T2a)
// ===========================================================================
//
// At a catalogue-good prime p (B = Algebra(O) unramified at p and O p-maximal,
// so O_p = M_2(Z_p)) every stable lattice a with p^e O <= a <= O is locally
// principal at p, and the stable lattices of index p^(2k) are in bijection with
// the index-p^k sublattices of Z_p^2, of which there are 1 + p + ... + p^k. We
// build Cat(p,e) by certified descent: from each node a take the p+1 maximal
// O-invariant subspaces of a/pa (each of codimension 2, since a/pa = S + S for
// the 2-dimensional simple M_2(F_p)-module S), pull them back, HNF-dedupe per
// depth, and ASSERT the deduped depth-k count is 1 + p + ... + p^k. The p-adic
// splitting is used only to CERTIFY goodness (IsUnramified/IspMaximal); the
// descent is exact F_p / Z linear algebra, so no p-adic precision enters the
// catalogue. Every depth-e entry is passed through the T2c filter; at a good
// prime it must be locally principal, so a filter failure is a hard error.
//
// Catalogue entry record: X (the 4x4 HNF O-coordinate matrix of the ideal a),
// n (= p^e, its reduced norm), Filter (the CMFilterFormat record from
// LocalPrincipalityData). Named distinctly from any T2b catalogue record so the
// parallel bad-prime section can coexist in this file; T2d may unify them.
CMGoodCatalogueEntry := recformat< X, n, Filter >;

// Deterministic total order on integer matrices: lexicographic on the row-major
// entries. Comparison function for Sort (returns -1, 0, or 1).
GoodCatMatLt := func< A, B |
    Eltseq(A) lt Eltseq(B) select -1 else (Eltseq(A) eq Eltseq(B) select 0 else 1) >;

// Preimage of an F_p-subspace W of a/pa as an HNF O-coordinate matrix. Wrows are
// integer lifts of an F_p-basis of W in the a-coordinates b_1..b_4 (rows of Xi);
// the preimage lattice in a-coordinates is lift(W) + p*Z^4 (index p^(4-dim W)),
// whose HNF basis is carried to O-coordinates by right multiplication by Xi.
function GoodCatPreimage(Wrows, Xi, p)
    d := #Wrows;
    lift := Matrix(Integers(), d, 4, &cat Wrows);
    stack := VerticalJoin(lift, p*IdentityMatrix(Integers(), 4));
    H := HermiteForm(stack);
    nz := [ H[l] : l in [1..Nrows(H)] | not IsZero(H[l]) ];
    assert #nz eq 4;
    return HermiteForm(Matrix(nz) * Xi);
end function;

// All 2-dimensional (codimension-2) subspaces of F_p^4, each as its unique 2x4
// reduced row-echelon matrix, in a fixed deterministic order. There are
// (p^2+1)(p^2+p+1) of them. Used only by the independent depth-1 cross-check.
function GoodCatAllTwoSubspaces(p)
    F := GF(p);
    res := [];
    for c1 in [1..4] do
      for c2 in [c1+1..4] do
        free1 := [ c : c in [c1+1..4] | c ne c2 ];   // row 1 non-pivot cols after c1
        free2 := [ c : c in [c2+1..4] ];              // row 2 non-pivot cols after c2
        nf := #free1 + #free2;
        for tup in CartesianPower([0..p-1], nf) do
            R := ZeroMatrix(F, 2, 4);
            R[1,c1] := 1; R[2,c2] := 1;
            j := 1;
            for c in free1 do R[1,c] := F ! tup[j]; j +:= 1; end for;
            for c in free2 do R[2,c] := F ! tup[j]; j +:= 1; end for;
            Append(~res, R);
        end for;
      end for;
    end for;
    return res;
end function;

// ---------------------------------------------------------------------------
// Exported: good-prime certification and catalogues.
// ---------------------------------------------------------------------------

intrinsic IsCatalogueGoodPrime(O::AlgQuatOrd, p::RngIntElt) -> BoolElt
{ True iff p is a catalogue-good prime for O: p is a rational prime, the algebra
  B = Algebra(O) is unramified at p, and O is p-maximal (IspMaximal(O, p)). Then
  O_p = M_2(Z_p) and GoodPrimeCatalogue applies. This is the catalogue notion of
  "good", broader than the theorem's good primes (odd, not dividing discrd(O)):
  p = 2 qualifies whenever B is unramified at 2 and O is 2-maximal. }
    require IsPrime(p): "p must be a rational prime";
    return IsUnramified(p, Algebra(O)) and IspMaximal(O, p);
end intrinsic;

intrinsic GoodPrimeCatalogue(O::AlgQuatOrd, p::RngIntElt, e::RngIntElt) -> SeqEnum
{ The catalogue Cat(p,e): every integral locally principal left O-ideal a with
  p^e O <= a <= O and [O:a] = p^(2e), at a catalogue-good prime p, built by
  certified descent through the p+1 maximal O-invariant subspaces of a/pa at each
  node (module machinery over the four left-multiplication matrices of Basis(O)),
  HNF-deduped per depth, with the count assert 1 + p + ... + p^k at every depth
  k <= e. Every depth-e entry passes the T2c filter (CMCandidatePremise at n = p^e,
  then LocalPrincipalityData with witnesses); a filter failure is a hard error,
  since at a good prime every such stable lattice is locally principal. Returns a
  deterministically ordered (lexicographic on HNF entries) SeqEnum of
  CMGoodCatalogueEntry records (fields X, n, Filter). Requires
  IsCatalogueGoodPrime(O, p); otherwise a clean error points at the bad-prime
  catalogue (T2b). }
    require IsPrime(p): "p must be a rational prime";
    require e ge 1: "e must be a positive integer";
    require IsCatalogueGoodPrime(O, p):
        "p is not a catalogue-good prime for O (B ramified at p, or O not p-maximal); "
        cat "use the bad-prime catalogue (T2b) for this prime";

    F := GF(p);
    I4 := IdentityMatrix(Integers(), 4);
    current := { I4 };   // depth 0: the unit ideal O
    for k in [1..e] do
        nextset := {};
        for X in current do
            Xi, _, AA := CMCandidatePremise(O, X, p^(k-1));   // A_i on this node a
            AAp := [ ChangeRing(AA[i], F) : i in [1..4] ];
            M := RModule(AAp);
            ms := MaximalSubmodules(M);
            // O_p = M_2(Z_p) forces a/pa = S + S: exactly p+1 maximal submodules,
            // each of codimension 2 (each descent step multiplies index by p^2).
            error if #ms ne p+1,
                "good-prime descent invariant violated at p=" cat IntegerToString(p)
                cat ": expected " cat IntegerToString(p+1)
                cat " maximal O-invariant subspaces of a/pa but found " cat IntegerToString(#ms);
            for S in ms do
                assert Dimension(S) eq 2;
                mor := Morphism(S, M);
                Wrows := [ [ Integers() ! x : x in Eltseq(mor[r]) ] : r in [1..Nrows(mor)] ];
                Xnew := GoodCatPreimage(Wrows, Xi, p);
                assert Determinant(Xnew) eq p^(2*k);   // stable lattice of index p^(2k)
                Include(~nextset, Xnew);
            end for;
        end for;
        current := nextset;
        expected := &+[ p^j : j in [0..k] ];   // 1 + p + ... + p^k
        error if #current ne expected,
            "good-prime catalogue count mismatch at p=" cat IntegerToString(p) cat ", depth "
            cat IntegerToString(k) cat ": expected " cat IntegerToString(expected)
            cat " deduped stable lattices but found " cat IntegerToString(#current);
    end for;

    norm := p^e;
    catalogue := [];
    for X in Sort([ M : M in current ], GoodCatMatLt) do
        Xi := CMCandidatePremise(O, X, norm);    // premise (index, containment, stability)
        assert Determinant(Xi) eq norm^2;        // norm-index identity [O:a] = n^2
        fr := LocalPrincipalityData(O, X, norm : Witnesses := true);
        error if not fr`LocallyPrincipal,
            "GOOD-PRIME FILTER FAILURE at p=" cat IntegerToString(p) cat ", e="
            cat IntegerToString(e) cat " (nu=" cat IntegerToString(fr`nu) cat ", n="
            cat IntegerToString(norm) cat "): a good-prime stable lattice must be locally "
            cat "principal, so this is a mathematical inconsistency";
        Append(~catalogue, rec< CMGoodCatalogueEntry | X := X, n := norm, Filter := fr >);
    end for;
    return catalogue;
end intrinsic;

intrinsic GoodPrimeDepth1CrossCheck(O::AlgQuatOrd, p::RngIntElt) -> BoolElt
{ Independent verification of Cat(p,1) at a catalogue-good prime p: enumerate ALL
  (p^2+1)(p^2+p+1) codimension-2 subspaces of O/pO by reduced row echelon form,
  keep those invariant under the four left-multiplication matrices of Basis(O),
  pull each back to an HNF O-coordinate lattice, and ASSERT the resulting set of
  HNF matrices equals that of GoodPrimeCatalogue(O, p, 1). This uses no
  maximal-submodule machinery, so it independently cross-checks the descent.
  Returns true on success; a mismatch is a hard error. }
    require IsPrime(p): "p must be a rational prime";
    require IsCatalogueGoodPrime(O, p):
        "p is not a catalogue-good prime for O; the depth-1 cross-check assumes O_p = M_2(Z_p)";
    F := GF(p);
    I4 := IdentityMatrix(Integers(), 4);
    _, _, AA := CMCandidatePremise(O, I4, 1);
    AAp := [ ChangeRing(AA[i], F) : i in [1..4] ];

    invset := {};
    for R in GoodCatAllTwoSubspaces(p) do
        good := true;
        for i in [1..4] do
            // W invariant under right multiplication by AAp[i] iff no rank gain.
            if Rank(VerticalJoin(R, R*AAp[i])) ne 2 then good := false; break; end if;
        end for;
        if good then
            Wrows := [ [ Integers() ! x : x in Eltseq(R[r]) ] : r in [1..2] ];
            Include(~invset, GoodCatPreimage(Wrows, I4, p));
        end if;
    end for;

    catset := { entry`X : entry in GoodPrimeCatalogue(O, p, 1) };
    error if invset ne catset,
        "depth-1 cross-check MISMATCH at p=" cat IntegerToString(p)
        cat ": independent invariant-subspace enumeration gave " cat IntegerToString(#invset)
        cat " lattices but the catalogue gave " cat IntegerToString(#catset);
    return true;
end intrinsic;

// ===========================================================================
// BAD-PRIME DESCENT CATALOGUES (T2b)
//
// For a prime p dividing discrd(O), the catalogue Cat(p, e) of locally principal
// candidate ideals of reduced norm p^e is built by bounded BFS descent through
// ALL maximal invariant subspaces of J/pJ starting from J = O, RETAINING every
// stable intermediate lattice (including non-locally-principal ones: pruning by
// local principality mid-descent is UNSOUND because a non-principal intermediate
// can have locally principal descendants). Only the terminal set (stable, index
// p^(2e), containing p^e O) is filtered through the T2c local-principality layer.
//
// Completeness is Nakayama: a proper O-stable sublattice K of a stable lattice J
// maps to a PROPER invariant subspace of J/pJ (else K + pJ = J forces K = J),
// hence lies in a MAXIMAL invariant subspace whose preimage is a stable
// intermediate above K; iterating reaches every stable lattice of bounded index.
// At bad p the algebra O/pO need not be semisimple, so maximal invariant
// subspaces have codimension 1 or 2 and the index steps by p or p^2:
// log_p[O:J'] is tracked exactly (never assumed even).
//
// Two prunes are sound and give termination: J' not containing p^e O is dropped
// (no sublattice of it contains p^e O either), and J' past index p^(2e) is
// dropped (every target lies between p^e O and O, so within (Z/p^e)^4). A node
// cap turns runaway growth into a clean error; a cutoff NEVER returns a partial
// catalogue.
//
// Reuses the T2c file-local helpers OrderMultData, ValidateStable, NuFromBasis,
// IdealBasisElts and the exported CMCandidatePremise / LocalPrincipalityData. The
// maximal-invariant-subspace helper is self-contained here (Bad prefix; the
// good-prime section keeps its own).
// ===========================================================================

// Catalogue entry: X (4x4 integral HNF, rows a Z-basis of the ideal in
// Basis(O)-coordinates), reduced norm n = p^e, and the T2c filter record.
BadCatalogueEntryFormat := recformat< X, n, Filter >;

// Descent statistics (audit trail plus the retained-intermediate evidence).
//   NodesExpanded     : nodes whose maximal invariant subspaces were computed.
//   NodesSeenPerIndex : SeqEnum of length 2e+1; entry k+1 is the count of distinct
//                       seen lattices with log_p[O:.] = k, for k = 0..2e.
//   TerminalCount     : pre-filter terminal candidates (index p^(2e), >= p^e O).
//   KeptCount         : locally principal terminals retained (= #entries).
//   KeptBelowNonLP    : kept entries with at least one non-locally-principal
//                       stable lattice strictly between O and them in the
//                       explored poset (evidence that retaining non-principal
//                       intermediates is necessary for completeness).
//   Codims            : sorted distinct codimensions encountered in descent.
//   Terminals         : the pre-filter terminal HNF matrices (for cross-checks).
BadCatalogueStatsFormat := recformat<
    p, e, NodesExpanded, NodesSeenPerIndex, TerminalCount, KeptCount,
    KeptBelowNonLP, Codims, Terminals >;

// Deterministic total order on integer matrices: lexicographic on Eltseq.
function BadHNFCmp(a, b)
    ea := Eltseq(a); eb := Eltseq(b);
    for i in [1..#ea] do
        if ea[i] lt eb[i] then return -1; end if;
        if ea[i] gt eb[i] then return 1; end if;
    end for;
    return 0;
end function;

// Maximal invariant subspaces of F_p^4 under RIGHT multiplication by the four
// matrices AAp (left multiplication by Basis(O) on the current node, mod p),
// returned as canonical row-echelon basis matrices, sorted. Canonicalizing and
// sorting makes the output independent of MeatAxe internal basis/order choices.
function BadMaxInvariantSubspaces(AAp)
    Mod := RModule(AAp);
    subs := MaximalSubmodules(Mod);
    bases := [ EchelonForm(Morphism(S, Mod)) : S in subs ];
    Sort(~bases, BadHNFCmp);
    return bases;
end function;

// Preimage in O-coordinates of a submodule W (row-echelon basis Wb over F_p in
// the J-basis) of J/pJ, where J has HNF matrix Xi (rows = Z-basis of J in
// Basis(O)-coordinates). The preimage is the Z-span of integer lifts of Wb and p
// times each J-basis vector, carried to O-coordinates by Xi, Hermite-reduced.
function BadPreimageLattice(Xi, Wb, p)
    d := Nrows(Wb);
    rowsJ := [ [ Integers()!Wb[a,b] : b in [1..4] ] : a in [1..d] ]
             cat [ [ (i eq j select p else 0) : j in [1..4] ] : i in [1..4] ];
    H := HermiteForm(Matrix(Integers(), rowsJ) * Xi);
    nz := [ H[l] : l in [1..Nrows(H)] | not IsZero(H[l]) ];
    return HermiteForm(Matrix(nz));
end function;

// Local principality at p of a stable lattice M (HNF Xi) known to contain p^e O
// (so [O:M] is a power of p). Principal iff that power is even, = p^(2k), and the
// content nu(M) = p^k. E is Basis(O) as algebra elements.
function BadLatLocPrincipal(E, Xi, p)
    g := Valuation(Determinant(Xi), p);
    if IsOdd(g) then return false; end if;
    return NuFromBasis(IdealBasisElts(E, Xi)) eq p^(g div 2);
end function;

intrinsic BadPrimeCatalogue(O::AlgQuatOrd, p::RngIntElt, e::RngIntElt) -> SeqEnum, Rec
{ Catalogue Cat(p, e) at a BAD prime p (dividing the reduced discriminant
  discrd(O)) and exponent e at least 1: all locally principal integral left
  O-ideals of reduced norm p^e, i.e. stable lattices J' with p^e O inside J'
  inside O, index [O:J'] = p^(2e), and J' principal over O_p. Built by bounded BFS
  from J = O through every maximal invariant subspace of J/pJ (right action of the
  four Basis(O) left-multiplication matrices mod p), retaining non-locally-
  principal intermediates, tracking log_p[O:J'] exactly, pruning any node that
  fails to contain p^e O or exceeds index p^(2e), and HNF-deduplicating. Terminal
  candidates pass through CMCandidatePremise and LocalPrincipalityData; only the
  locally principal ones are kept (rejection is expected at bad primes). Requires
  p prime dividing discrd(O); a good prime redirects to the good-prime catalogue.
  Returns the kept entries (BadCatalogueEntryFormat records X, n = p^e, filter
  record; ordered canonically by HNF) and a BadCatalogueStatsFormat statistics
  record. Raises a clean error on a node-cap overflow rather than returning a
  partial catalogue. }
    require IsPrime(p): "p must be a prime";
    require e ge 1: "e must be a positive integer";
    require IsDivisibleBy(Integers()!Discriminant(O), p):
        "p must divide discrd(O); a good prime (p unramified in the algebra, O p-maximal) uses the good-prime catalogue (T2a)";

    _, E, _, _, AO := OrderMultData(O);
    F := GF(p);
    pe := p^e;
    twoe := 2*e;
    // One matrix parent for every lattice, so set/AssociativeArray keys compare by
    // value (IdentityMatrix and HermiteForm can otherwise land in distinct parents).
    MS := RMatrixSpace(Integers(), 4, 4);
    I4 := MS ! IdentityMatrix(Integers(), 4);
    NODECAP := 200000;

    // BFS state. seen dedupes; logOf carries log_p[O:.]; frontier holds nodes
    // still to expand (log < 2e); terminals holds leaves (log = 2e, >= p^e O).
    seen := { I4 };
    logOf := AssociativeArray();
    logOf[I4] := 0;
    seenList := [ I4 ];
    frontier := [ I4 ];
    terminals := [ MS | ];
    nodesExpanded := 0;
    codimSet := { Integers() | };

    while #frontier gt 0 do
        Xi := frontier[1];
        Remove(~frontier, 1);
        g := logOf[Xi];
        nodesExpanded +:= 1;

        ok, _, AA := ValidateStable(Xi, AO);
        assert ok;                          // every node is O-stable by construction
        AAp := [ ChangeRing(AA[i], F) : i in [1..4] ];

        newnodes := [ MS | ];
        for Wb in BadMaxInvariantSubspaces(AAp) do
            codim := 4 - Nrows(Wb);
            Include(~codimSet, codim);
            Xp := MS ! BadPreimageLattice(Xi, Wb, p);
            gp := g + codim;
            assert Determinant(Xp) eq p^gp;         // exact index cross-check
            if gp gt twoe then continue; end if;    // stop past index p^(2e)
            XpQ := ChangeRing(Xp, Rationals());
            if not &and[Denominator(c) eq 1 : c in Eltseq(pe*XpQ^-1)] then
                continue;                            // does not contain p^e O: sound prune
            end if;
            if Xp in seen then continue; end if;
            Include(~seen, Xp);
            logOf[Xp] := gp;
            Append(~seenList, Xp);
            Append(~newnodes, Xp);
            require #seen le NODECAP:
                "BadPrimeCatalogue: node cap exceeded; refusing to return a partial catalogue";
        end for;

        Sort(~newnodes, BadHNFCmp);
        for Xp in newnodes do
            if logOf[Xp] eq twoe then
                Append(~terminals, Xp);
            else
                Append(~frontier, Xp);
            end if;
        end for;
    end while;

    Sort(~terminals, BadHNFCmp);

    // Filter terminals through the T2c local-principality layer; keep principal.
    entries := [ ];
    for Xp in terminals do
        _ := CMCandidatePremise(O, Xp, pe);
        fdata := LocalPrincipalityData(O, Xp, pe : Witnesses := true);
        if fdata`LocallyPrincipal then
            Append(~entries, rec< BadCatalogueEntryFormat |
                X := Xp, n := pe, Filter := fdata >);
        end if;
    end for;

    // Evidence: kept entries sitting below a non-locally-principal intermediate.
    keptbelow := 0;
    for entry in entries do
        Xk := entry`X;
        XkQ := ChangeRing(Xk, Rationals());
        below := false;
        for M in seenList do
            if M eq Xk or M eq I4 then continue; end if;      // strict, and not O
            T := XkQ * (ChangeRing(M, Rationals()))^-1;        // Xk <= M iff integral
            if not &and[Denominator(c) eq 1 : c in Eltseq(T)] then continue; end if;
            if not BadLatLocPrincipal(E, M, p) then below := true; break; end if;
        end for;
        if below then keptbelow +:= 1; end if;
    end for;

    perindex := [ 0 : k in [0..twoe] ];
    for M in seenList do
        perindex[logOf[M] + 1] +:= 1;
    end for;

    stats := rec< BadCatalogueStatsFormat |
        p := p, e := e,
        NodesExpanded := nodesExpanded,
        NodesSeenPerIndex := perindex,
        TerminalCount := #terminals,
        KeptCount := #entries,
        KeptBelowNonLP := keptbelow,
        Codims := Sort(SetToSequence(codimSet)),
        Terminals := terminals >;

    return entries, stats;
end intrinsic;
