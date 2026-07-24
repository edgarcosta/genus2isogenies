/*
 * Quadratic-lattice layer for the CM gluing-levels engine (magma/cmlevels).
 *
 * Given two full-rank lattices I, J in a common imaginary quadratic field L,
 * this file computes, by exact rational 2x2 linear algebra over the fixed
 * Q-basis [1, L.1] of L (no floating point, no non-maximal-order ideal
 * arithmetic), the objects the order layer consumes:
 *
 *   R  = (I:I) = {x in L : x*I subset I}   (multiplicator ring of I)
 *   R' = (J:J)
 *   M  = (I:J) = {x in L : x*J subset I}   (= Hom(E', E))
 *   R0 = R meet R'                          (an order, Z-basis [1, r])
 *   delta = covol(J)/covol(I) in Q_{>0}     (= |det B_J| / |det B_I|)
 *
 * A lattice is represented internally by a 2x2 rational matrix whose rows are
 * the coordinate vectors (via Eltseq, i.e. the power basis [1, L.1]) of a
 * Z-basis; publicly it is a length-2 SeqEnum of field elements. Every derived
 * basis is put in Hermite normal form so it is a function of the lattice, not
 * of the presentation supplied by the caller. The common denominator of a
 * basis matrix equals min{d : d*Lat subset Z^2}, an invariant of the lattice,
 * so the HNF below is well defined on any generating set of a fixed lattice.
 */

// ---------------------------------------------------------------------------
// Internal helpers (file-local; matrix = rational row-basis of a lattice).
// ---------------------------------------------------------------------------

// Coordinate row vector of a field element in the power basis [1, L.1].
function CoordVec(x)
    return Vector(Rationals(), Eltseq(x));
end function;

// 2x2 rational basis matrix (rows = coordinates of the two basis elements).
function BasisMat(bas)
    return Matrix(Rationals(), [Eltseq(bas[1]), Eltseq(bas[2])]);
end function;

// Reconstruct the field elements of the rows of a 2x2 matrix P.
function MatToElts(P, L)
    return [L ! Eltseq(P[1]), L ! Eltseq(P[2])];
end function;

// Is the rational row vector v in the Z-row-span of the rational basis B?
function InLat(v, B)
    w := v * B^-1;
    return &and[Denominator(c) eq 1 : c in Eltseq(w)];
end function;

// Canonical HNF row-basis of the Z-lattice spanned by the rows of a rational
// matrix M (any number of rows, rank 2). Returns a 2x2 rational matrix.
function LatHNF(M)
    d := Denominator(M);
    H := HermiteForm(Matrix(Integers(), d*M));
    B := ChangeRing(Submatrix(H, 1, 1, 2, 2), Rationals()) / d;
    assert Determinant(B) ne 0;
    return B;
end function;

// Dual lattice w.r.t. the standard pairing: {x : <x, row_i(P)> in Z for all i}.
function LatDual(P)
    return LatHNF(Transpose(P^-1));
end function;

// Sum (Z-span of the union of two bases) and intersection of two lattices.
function LatSum(P, Q)
    return LatHNF(VerticalJoin(P, Q));
end function;

function LatMeet(P, Q)
    return LatDual(LatSum(LatDual(P), LatDual(Q)));
end function;

// Colon lattice (I:J) = {x in L : x*J subset I}, from row-basis matrices of I, J.
// Condition on coords x: for each Z-basis element jb of J, x*RepMat(jb) lands in
// the Z-row-span of B_I (verified convention: coords(x*jb) = coords(x)*RepMat(jb)).
// Stacking both columns, x*C in Z^4, so (I:J) is the dual of the column lattice of C.
function LatColon(BI, BJ, L)
    BIinv := BI^-1;
    j1 := L ! Eltseq(BJ[1]);
    j2 := L ! Eltseq(BJ[2]);
    R1 := ChangeRing(RepresentationMatrix(j1), Rationals());
    R2 := ChangeRing(RepresentationMatrix(j2), Rationals());
    C := HorizontalJoin(R1 * BIinv, R2 * BIinv);
    return LatDual(LatHNF(Transpose(C)));
end function;

// Normalize a caller-supplied lattice basis to field elements of L, returning
// L and the length-2 list of coerced elements. Raises clean errors on bad input.
function NormalizeBasis(bas, name)
    error if #bas ne 2, name * " must be a 2-element basis (got " * IntegerToString(#bas) * " elements)";
    U := Universe(bas);
    if ISA(Type(U), RngOrd) then
        L := NumberField(U);
    else
        L := U;
    end if;
    error if not ISA(Type(L), FldNum), name * " must consist of number field elements";
    return L, [L ! x : x in bas];
end function;

// Validate a pair (I, J) as full-rank lattices in a common imaginary quadratic
// field. Returns L, the coerced bases, and their rational basis matrices.
function ValidatePair(Ibas, Jbas)
    LI, Ie := NormalizeBasis(Ibas, "I");
    LJ, Je := NormalizeBasis(Jbas, "J");
    error if LI ne LJ, "I and J must lie in the same number field";
    L := LI;
    error if Degree(L) ne 2, "L must be a quadratic field (degree 2)";
    s1, s2 := Signature(L);
    error if s1 ne 0, "L must be imaginary quadratic (found a real embedding)";
    BI := BasisMat(Ie);
    BJ := BasisMat(Je);
    error if Determinant(BI) eq 0, "I is rank deficient (basis elements are Q-linearly dependent)";
    error if Determinant(BJ) eq 0, "J is rank deficient (basis elements are Q-linearly dependent)";
    return L, Ie, Je, BI, BJ;
end function;

// ---------------------------------------------------------------------------
// Exported intrinsics.
// ---------------------------------------------------------------------------

intrinsic CanonicalLatticeBasis(bas::SeqEnum) -> SeqEnum
{ The canonical (Hermite normal form) Z-basis of the full-rank lattice in an
  imaginary quadratic field spanned by bas; a function of the lattice, not of
  the presentation. }
    L, e := NormalizeBasis(bas, "basis");
    B := BasisMat(e);
    require Determinant(B) ne 0: "basis is rank deficient";
    return MatToElts(LatHNF(B), L);
end intrinsic;

intrinsic LatticeColon(Ibas::SeqEnum, Jbas::SeqEnum) -> SeqEnum
{ Canonical Z-basis of the colon lattice (I:J), the x in L with x*J contained in
  I, for full-rank lattices I, J given by Z-bases in a common imaginary quadratic
  field. }
    L, _, _, BI, BJ := ValidatePair(Ibas, Jbas);
    return MatToElts(LatColon(BI, BJ, L), L);
end intrinsic;

intrinsic LatticeMeet(Ibas::SeqEnum, Jbas::SeqEnum) -> SeqEnum
{ Canonical Z-basis of the intersection lattice I meet J, for full-rank lattices
  I, J given by Z-bases in a common imaginary quadratic field. }
    L, _, _, BI, BJ := ValidatePair(Ibas, Jbas);
    return MatToElts(LatMeet(BI, BJ), L);
end intrinsic;

intrinsic LatticeCovolumeRatio(Ibas::SeqEnum, Jbas::SeqEnum) -> FldRatElt
{ The exact positive-rational covolume ratio delta = covol(J)/covol(I) =
  |det B_J|/|det B_I|, computed algebraically over the fixed Q-basis of L. }
    _, _, _, BI, BJ := ValidatePair(Ibas, Jbas);
    return AbsoluteValue(Determinant(BJ)) / AbsoluteValue(Determinant(BI));
end intrinsic;

intrinsic CMLatticeConstruction(Ibas::SeqEnum, Jbas::SeqEnum)
    -> SeqEnum, SeqEnum, SeqEnum, SeqEnum, FldRatElt
{ Given Z-bases of full-rank lattices I, J in a common imaginary quadratic field
  L, return canonical Z-bases of R = (I:I), R' = (J:J), R0 = R meet R', M = (I:J),
  together with the positive rational delta = covol(J)/covol(I). All lattice-layer
  invariants are asserted before returning. }
    L, Ie, Je, BI, BJ := ValidatePair(Ibas, Jbas);

    Rm  := LatColon(BI, BI, L);
    Rpm := LatColon(BJ, BJ, L);
    Mm  := LatColon(BI, BJ, L);
    R0m := LatMeet(Rm, Rpm);
    delta := AbsoluteValue(Determinant(BJ)) / AbsoluteValue(Determinant(BI));

    Rel   := MatToElts(Rm, L);
    Rpel  := MatToElts(Rpm, L);
    R0el  := MatToElts(R0m, L);
    Mel   := MatToElts(Mm, L);
    one := L ! 1;

    // 1 lies in each multiplicator ring / order.
    assert InLat(CoordVec(one), Rm);
    assert InLat(CoordVec(one), Rpm);
    assert InLat(CoordVec(one), R0m);

    // Each of R, R', R0 is closed under multiplication.
    assert forall{ <a, b> : a in Rel,  b in Rel  | InLat(CoordVec(a*b), Rm)  };
    assert forall{ <a, b> : a in Rpel, b in Rpel | InLat(CoordVec(a*b), Rpm) };
    assert forall{ <a, b> : a in R0el, b in R0el | InLat(CoordVec(a*b), R0m) };

    // Module properties: R*I subset I, R'*J subset J, M*J subset I, R0*M subset M.
    assert forall{ <r, x> : r in Rel,  x in Ie | InLat(CoordVec(r*x), BI) };
    assert forall{ <r, y> : r in Rpel, y in Je | InLat(CoordVec(r*y), BJ) };
    assert forall{ <m, y> : m in Mel,  y in Je | InLat(CoordVec(m*y), BI) };
    assert forall{ <u, m> : u in R0el, m in Mel | InLat(CoordVec(u*m), Mm) };

    // Orientation certificate for delta: (J:I) = delta*conj(M) as lattices,
    // with M = (I:J). An inverted delta breaks this, since the two sides
    // scale oppositely in delta.
    Mconj := LatHNF(BasisMat([Trace(m) - m : m in Mel]));
    assert LatColon(BJ, BI, L) eq LatHNF(delta * Mconj);

    assert delta gt 0;

    return Rel, Rpel, R0el, Mel, delta;
end intrinsic;
