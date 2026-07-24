/*
 * Equivalence layer for the CM gluing-levels engine (magma/cmlevels).
 *
 * Groups a stream of integral, locally principal left O-ideals (O the definite
 * quaternion order from CMSquareOrder, task T1) into left-ideal classes with
 * certified minimal representatives. The whole layer rests on ONE SVP
 * normalization (ReducedNormMinimum): for a Z-basis c_1..c_4 of a lattice C in
 * B = Algebra(O), the polar Gram G_ij = Trd(c_i conj(c_j)) satisfies
 * x G x^t = 2 nrd(sum x_i c_i); clearing denominators by d and calling
 * Minimum(LatticeWithGram(d*G) : Proof := true) certifies 2*d*(min nrd on C).
 * Every reduced norm below divides out the 2*d factor explicitly.
 *
 * Equivalence convention: a ~ b iff b = a*gamma for gamma in B^* (right
 * multiplication; both are left O-ideals with left order O). The certified test
 * forms C = conj(a)*b (16 products, HNF to rank 4): a ~ b iff C contains z with
 * nrd(z) = nrd(a)*nrd(b), and then gamma = z/nrd(a) with a*gamma = b re-verified
 * by exact HNF lattice equality. Stock Magma only; exact arithmetic throughout.
 *
 * Ideal representation (shared with enumerate.m, by orchestrator agreement): an
 * integral left O-ideal a of reduced norm n is its 4x4 integral HNF coordinate
 * matrix X w.r.t. Basis(O) (rows = Z-basis of a in O-coordinates), plus n, with
 * det(X) = n^2. This file does NOT depend on enumerate.m; premise checks are
 * inline. Local principality is an input CONTRACT (T2c/T2d filter for it): the
 * certified minimum-equals-product characterization needs it, and the
 * min-nrd >= n_a*n_b assert below fails loudly if a non-invertible ideal slips in.
 */

// Class record produced by ClassifyByEquivalence.
//   Order      : the order O (so AssertFirstOccurrenceIsGlobalMin(cls) is self-contained)
//   RepX, RepN : minimal integral representative (first occurrence) and its reduced norm
//   RepIndex   : index of the representative in the input list
//   Members    : indices (into the input list) of all ideals in this class
//   Gammas     : parallel to Members; Gammas[t] is the exact, re-verified gamma with
//                member = Rep * Gammas[t] (Gammas at the representative is 1)
//   BucketKey  : the exact invariant key the representative bucketed under
CMIdealClass := recformat< Order, RepX, RepN, RepIndex, Members, Gammas, BucketKey >;

// Replay certificate for GlobalMinimalNorm.
//   NMin            : the certified global minimal integral norm of the class
//   ScaleD          : denominator-clearing factor d of the SVP lattice
//   LatticeMinimum  : m = Minimum(LatticeWithGram(d*G)) = 2*d*NMin/n_a
//   ShortestVectors : certified shortest vectors of conj(a)*O/n_a, as elements of B
CMGlobalMinCert := recformat< NMin, ScaleD, LatticeMinimum, ShortestVectors >;

// ---------------------------------------------------------------------------
// File-local helpers.
// ---------------------------------------------------------------------------

// True iff every entry of the rational matrix M is a rational integer.
function IsIntegralMat(M)
    return &and[ Denominator(c) eq 1 : c in Eltseq(M) ];
end function;

// Coordinate context of O: the algebra B, the O-basis as B-elements (qb), the
// row-coordinate matrix coO (rows = Eltseq of the basis) and its inverse, and the
// four integral left-multiplication matrices Lmat[i] with
// Ocoords(o_i * x) = Ocoords(x) * Lmat[i]. Integrality of Lmat is exactly closure
// of O under multiplication and is asserted here (a sanity check on O).
function OrderData(O)
    B := Algebra(O);
    qb := [ B ! b : b in Basis(O) ];
    coO := Matrix(Rationals(), [ Eltseq(x) : x in qb ]);
    coOinv := coO^-1;
    Lmat := [];
    for ii in [1..4] do
        rows := [ Vector(Rationals(), Eltseq(qb[ii]*qb[kk]))*coOinv : kk in [1..4] ];
        M := Matrix(Rationals(), [ Eltseq(r) : r in rows ]);
        assert IsIntegralMat(M);
        Append(~Lmat, M);
    end for;
    return B, qb, coO, coOinv, Lmat;
end function;

// Build the B-elements of a Z-basis from a rational O-coordinate matrix Xr.
function RowsToElements(qb, Xr)
    return [ &+[ Xr[k][l]*qb[l] : l in [1..4] ] : k in [1..4] ];
end function;

// Shared premise check. Errors (catchable) on any violation; on success returns
// B, qb, coO, coOinv, Lmat and the rational form Xr of X.
function CheckPremise(O, X, n)
    error if Nrows(X) ne 4 or Ncols(X) ne 4, "ideal matrix X must be 4x4";
    error if n le 0, "claimed reduced norm n must be a positive integer";
    Xr := Matrix(Rationals(), X);
    error if not IsIntegralMat(Xr),
        "ideal matrix X must be integral (rows are O-coordinates of a Z-basis of a <= O)";
    dt := Determinant(Xr);
    error if dt eq 0, "ideal matrix X must be full rank";
    error if dt ne n^2, Sprintf("det(X) = %o must equal n^2 = %o (norm-index premise [O:a] = n^2)", dt, n^2);
    B, qb, coO, coOinv, Lmat := OrderData(O);
    Xinv := Xr^-1;
    error if not IsIntegralMat(n*Xinv), "n*O must be contained in a (n*X^-1 is not integral)";
    for ii in [1..4] do
        error if not IsIntegralMat(Xr*Lmat[ii]*Xinv),
            Sprintf("row lattice is not left O-stable (fails at O-basis element %o)", ii);
    end for;
    return B, qb, coO, coOinv, Lmat, Xr;
end function;

// HNF O-coordinate matrix (4x4 integral) of the rank-4 lattice spanned by the
// B-elements gens; errors if they do not lie in O or do not span rank 4.
function CoordMatrix(qb, coOinv, gens)
    B := Universe(qb);
    rows := [ Vector(Rationals(), Eltseq(B!x))*coOinv : x in gens ];
    M := Matrix(Rationals(), [ Eltseq(r) : r in rows ]);
    error if not IsIntegralMat(M), "generators must lie in O (integral O-coordinates)";
    H := HermiteForm(Matrix(Integers(), M));
    nz := [ H[t] : t in [1..Nrows(H)] | not IsZero(H[t]) ];
    error if #nz ne 4, Sprintf("generators must span a rank-4 lattice (got rank %o)", #nz);
    return Matrix(Integers(), [ Eltseq(r) : r in nz ]);
end function;

// ---------------------------------------------------------------------------
// Exported intrinsics: converters and premise.
// ---------------------------------------------------------------------------

intrinsic IdealBasisElements(O::AlgQuatOrd, X::Mtrx) -> SeqEnum
{ The Z-basis of the left O-ideal with O-coordinate matrix X (rows = O-coordinates),
  as a sequence of elements of Algebra(O). }
    require Nrows(X) eq 4 and Ncols(X) eq 4: "X must be 4x4";
    _, qb := OrderData(O);
    return RowsToElements(qb, Matrix(Rationals(), X));
end intrinsic;

intrinsic IdealCoordinateMatrix(O::AlgQuatOrd, gens::SeqEnum) -> Mtrx
{ The 4x4 integral HNF O-coordinate matrix of the rank-4 lattice spanned (over Z)
  by gens, a sequence of elements of O (or coercible into Algebra(O)). }
    _, qb, _, coOinv := OrderData(O);
    return CoordMatrix(qb, coOinv, gens);
end intrinsic;

intrinsic AssertIdealPremise(O::AlgQuatOrd, X::Mtrx, n::RngIntElt) -> BoolElt
{ Assert the ideal premise for (X, n): X integral 4x4, full rank, det(X) = n^2,
  n*O contained in the row lattice, and left O-stability. Errors on any violation;
  returns true on success. }
    _ := CheckPremise(O, X, n);
    return true;
end intrinsic;

// ---------------------------------------------------------------------------
// Exported intrinsics: the SVP normalization and the product lattice.
// ---------------------------------------------------------------------------

intrinsic ReducedNormMinimum(C::SeqEnum) -> FldRatElt, SeqEnum, RngIntElt, RngIntElt
{ The single SVP normalization helper. Given a Z-basis C = [c_1..c_4] of a rank-4
  lattice in a definite quaternion algebra B, form the polar Gram
  G_ij = Trd(c_i conj(c_j)); with d clearing its denominators and
  L = LatticeWithGram(d*G), m = Minimum(L : Proof := true) certifies m = 2*d*(min
  nrd on C). Returns (min nrd on C) = m/(2*d) as an exact rational, the certified
  shortest vectors as elements of B, and (d, m) for replay. }
    require #C eq 4: "C must be a Z-basis of a rank-4 lattice (exactly 4 elements)";
    G := Matrix(Rationals(), 4, 4,
        [[ Trace(C[a]*Conjugate(C[b])) : b in [1..4] ] : a in [1..4] ]);
    require G eq Transpose(G): "polar Gram must be symmetric";
    d := LCM([ Denominator(c) : c in Eltseq(G) ] cat [1]);
    dG := Matrix(Integers(), d*G);
    require IsPositiveDefinite(dG):
        "polar Gram must be positive definite (definite algebra, full-rank lattice)";
    L := LatticeWithGram(dG);
    m := Minimum(L : Proof := true);
    sv := ShortestVectors(L);
    zs := [ &+[ (Integers()!Eltseq(s)[l])*C[l] : l in [1..4] ] : s in sv ];
    minnrd := (Rationals()!m) / (2*d);
    return minnrd, zs, d, Integers()!m;
end intrinsic;

intrinsic ConjProductLattice(ab::SeqEnum, bb::SeqEnum) -> SeqEnum
{ A Z-basis of the rank-4 lattice conj(a)*b in B = Universe(ab): the Z-span of the
  16 products conj(x_i)*y_j for the given Z-bases ab of a and bb of b, reduced to
  a rank-4 basis by HNF. Errors if the span is not rank 4. }
    require #ab ge 1 and #bb ge 1: "both bases must be nonempty";
    B := Universe(ab);
    prods := [ Conjugate(ab[p])*bb[q] : p in [1..#ab], q in [1..#bb] ];
    P := Matrix(Rationals(), [ Eltseq(x) : x in prods ]);
    den := LCM([ Denominator(c) : c in Eltseq(P) ] cat [1]);
    H := HermiteForm(Matrix(Integers(), den*P));
    nz := [ H[t] : t in [1..Nrows(H)] | not IsZero(H[t]) ];
    require #nz eq 4: Sprintf("conj(a)*b must have rank 4 (got rank %o)", #nz);
    return [ B ! Eltseq(Vector(Rationals(), r)/den) : r in nz ];
end intrinsic;

// ---------------------------------------------------------------------------
// Exported intrinsics: equivalence, buckets, classes, minima.
// ---------------------------------------------------------------------------

intrinsic AreLeftEquivalent(O::AlgQuatOrd, Xa::Mtrx, na::RngIntElt, Xb::Mtrx, nb::RngIntElt)
    -> BoolElt, AlgQuatElt
{ Certified test of a ~ b for left O-ideals a = (Xa, na), b = (Xb, nb). Forms
  C = conj(a)*b, runs the SVP normalization, asserts min-nrd >= na*nb, and returns
  true iff min-nrd = na*nb. On equality returns gamma = z/na (z a certified shortest
  vector) after re-verifying a*gamma = b by exact HNF lattice equality; on
  inequivalence returns false and 0. The certified path is the SVP one. }
    B, qb, _, coOinv, _, Xar := CheckPremise(O, Xa, na);
    _, _, _, _, _, Xbr := CheckPremise(O, Xb, nb);
    abasis := RowsToElements(qb, Xar);
    bbasis := RowsToElements(qb, Xbr);
    C := ConjProductLattice(abasis, bbasis);
    minnrd, zs := ReducedNormMinimum(C);
    target := na*nb;
    error if minnrd lt target,
        Sprintf("SVP minimum %o is below na*nb = %o: a normalization bug, or a non-locally-principal ideal was supplied (contract violation)", minnrd, target);
    error if not IsIntegral(minnrd/target),
        Sprintf("SVP minimum %o is not a multiple of na*nb = %o: a normalization bug, or a non-locally-principal ideal was supplied", minnrd, target);
    if minnrd ne target then
        return false, B ! 0;
    end if;
    gamma := zs[1] / na;
    require Norm(gamma) eq nb/na: "recovered gamma has the wrong reduced norm (normalization bug)";
    Xag := CoordMatrix(qb, coOinv, [ x*gamma : x in abasis ]);
    require Xag eq HermiteForm(Matrix(Integers(), Xbr)):
        "a*gamma does not equal b on exact HNF re-verification (equivalence witness bug)";
    return true, gamma;
end intrinsic;

intrinsic BucketKey(O::AlgQuatOrd, X::Mtrx, n::RngIntElt) -> SeqEnum
{ A hashable exact class invariant for bucketing: Eltseq(ReducedGramMatrix(a))/n,
  a = lideal<O | Z-basis of X>. Equivalent ideals get equal keys (buckets only
  accelerate; inequivalent ideals may collide and are separated by AreLeftEquivalent).
  The lideal container is used only as a container: its basis is independently
  asserted to match X, and its reduced norm to match n. }
    B, qb, _, _, _, Xr := CheckPremise(O, X, n);
    basis := RowsToElements(qb, Xr);
    I := lideal< O | basis >;
    require Norm(I) eq n: "lideal container reduced norm disagrees with n";
    require IdealCoordinateMatrix(O, [ B ! x : x in Basis(I) ]) eq HermiteForm(Matrix(Integers(), Xr)):
        "lideal container basis disagrees with X";
    R := ReducedGramMatrix(I);
    return [ (Rationals()!R[a][b]) / n : a in [1..4], b in [1..4] ];
end intrinsic;

intrinsic ClassifyByEquivalence(O::AlgQuatOrd, ideals::SeqEnum) -> SeqEnum
{ Group ideals (a sequence of <X, n> tuples given in nondecreasing n) into left-ideal
  classes. Each incoming ideal is compared (by AreLeftEquivalent) only against the
  existing class representatives in its bucket; the first occurrence of a class is its
  minimal integral representative. Returns a sequence of CMIdealClass records. A final
  certification pass asserts that no two class representatives are equivalent, which
  certifies the bucketed grouping was complete (equivalent ideals never split across
  classes) independently of any canonicity assumption on the bucket key. }
    require Type(ideals) eq SeqEnum: "ideals must be a sequence of <X, n> tuples";
    for t in [2..#ideals] do
        require ideals[t][2] ge ideals[t-1][2]:
            "ideals must be given in nondecreasing reduced norm";
    end for;
    B := Algebra(O);
    classes := [];        // SeqEnum of CMIdealClass records
    bkeys := [];          // SeqEnum of bucket keys
    bclass := [];         // SeqEnum of SeqEnum[RngIntElt]: class indices per bucket
    for idx in [1..#ideals] do
        X := ideals[idx][1]; n := ideals[idx][2];
        key := BucketKey(O, X, n);
        bi := 0;
        for t in [1..#bkeys] do if bkeys[t] eq key then bi := t; break; end if; end for;
        found := false;
        if bi ne 0 then
            for ci in bclass[bi] do
                eqv, gamma := AreLeftEquivalent(O, classes[ci]`RepX, classes[ci]`RepN, X, n);
                if eqv then
                    c := classes[ci];
                    c`Members := Append(c`Members, idx);
                    c`Gammas := Append(c`Gammas, gamma);
                    classes[ci] := c;
                    found := true;
                    break;
                end if;
            end for;
        end if;
        if not found then
            newc := rec< CMIdealClass |
                Order := O, RepX := Matrix(Integers(), Matrix(Rationals(), X)),
                RepN := n, RepIndex := idx, Members := [idx], Gammas := [B ! 1],
                BucketKey := key >;
            Append(~classes, newc);
            if bi eq 0 then
                Append(~bkeys, key); Append(~bclass, [#classes]);
            else
                Append(~bclass[bi], #classes);
            end if;
        end if;
    end for;
    for a in [1..#classes] do
        for b in [a+1..#classes] do
            eqab := AreLeftEquivalent(O, classes[a]`RepX, classes[a]`RepN,
                                         classes[b]`RepX, classes[b]`RepN);
            error if eqab,
                "two class representatives are equivalent: the bucket key failed to be a class invariant and a class was split (classification incomplete); investigate ReducedGramMatrix canonicity";
        end for;
    end for;
    return classes;
end intrinsic;

intrinsic GlobalMinimalNorm(O::AlgQuatOrd, Xa::Mtrx, na::RngIntElt) -> RngIntElt, Rec
{ The global minimal integral norm of the left-ideal class of a = (Xa, na),
  independent of enumeration order: na * (min nrd on conj(a)*O / na). Returns the
  minimum and a CMGlobalMinCert replay certificate. }
    _, qb, _, _, _, Xar := CheckPremise(O, Xa, na);
    abasis := RowsToElements(qb, Xar);
    H := ConjProductLattice(abasis, qb);          // conj(a)*O
    Hs := [ h / na : h in H ];                     // scale by 1/na
    minH, zs, d, m := ReducedNormMinimum(Hs);
    nmin := na * minH;
    require IsIntegral(nmin) and nmin ge 1:
        Sprintf("global minimal norm %o must be a positive integer (normalization bug)", nmin);
    nmin := Integers() ! nmin;
    cert := rec< CMGlobalMinCert |
        NMin := nmin, ScaleD := d, LatticeMinimum := m, ShortestVectors := zs >;
    return nmin, cert;
end intrinsic;

intrinsic IsFirstOccurrenceGlobalMin(O::AlgQuatOrd, X::Mtrx, n::RngIntElt) -> BoolElt, RngIntElt
{ Boolean form of the first-occurrence guard: returns whether the global minimal norm
  of the class of (X, n) equals n, together with that global minimum. Use this to test
  the guard without raising. }
    nmin := GlobalMinimalNorm(O, X, n);
    return nmin eq n, nmin;
end intrinsic;

intrinsic AssertFirstOccurrenceIsGlobalMin(O::AlgQuatOrd, X::Mtrx, n::RngIntElt)
{ Hard-assert that n equals the global minimal norm of the class of (X, n). A failure
  means a smaller-norm equivalent ideal was missed, i.e. the enumeration feeding this
  representative was incomplete. }
    ok, nmin := IsFirstOccurrenceGlobalMin(O, X, n);
    error if not ok,
        Sprintf("first-occurrence norm %o exceeds the global class minimum %o: enumeration incomplete (a smaller-norm equivalent integral ideal was missed)", n, nmin);
end intrinsic;

intrinsic AssertFirstOccurrenceIsGlobalMin(cls::Rec)
{ As above, for a CMIdealClass record produced by ClassifyByEquivalence. }
    AssertFirstOccurrenceIsGlobalMin(cls`Order, cls`RepX, cls`RepN);
end intrinsic;
