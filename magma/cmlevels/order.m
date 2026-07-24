/*
 * Order-construction layer for the CM gluing-levels engine (magma/cmlevels).
 *
 * From the lattice data (R0 = R meet R', M = (I:J), delta) produced by
 * lattices.m, build the definite quaternion order
 *
 *     O = (R meet R') + M*j   inside   B = (s, -delta)_Q,
 *
 * where s = D_L < 0 is the discriminant of the imaginary quadratic field L and
 * tau = sqrt(D_L) is the deterministic trace-zero generator of L (documented
 * choice; the field discriminant is a canonical invariant of L, so O's algebra
 * presentation does not depend on how the caller wrote L). L embeds into B by
 * x + y*tau -> x + y*i (i^2 = s); Magma's convention j*i = -i*j then makes
 * j*a = conj(a)*j for a in the image of L, so nrd(a + b*j) = N_L(a) + delta*N_L(b).
 *
 * The order is built with the EXACT Z-basis [1, r, m1*j, m2*j] ([1, r] the
 * canonical basis of R0, [m1, m2] of M) via QuaternionOrder(Integers(), basis),
 * never the basis-only constructor which may replace the basis. Discriminant(O)
 * in Magma is the reduced discriminant. Four Layer-1 asserts run on every
 * construction (order closure, the nrd block formula, det(polar Gram) = discrd^2,
 * integrality of reduced trace/norm on the basis).
 */

// Record schema for the audit reason sets (section 2 decision 5): quantities
// whose odd prime divisors are provably subsumed by discrd(O); reported, not
// promoted to extra primes.
//   discR, discRp        : discriminants of R = (I:I), R' = (J:J)
//   DegreeFormContent    : content of the degree form delta*N_L restricted to M,
//                          gcd(q(m1), q(m2), q(m1+m2)-q(m1)-q(m2)), q(b)=delta*N_L(b)
//                          (cross-term form: no spurious factor 2 at p = 2)
//   OddPrimesSubsumed    : computed check that every odd prime dividing discR,
//                          discRp, or DegreeFormContent divides discrd(O)
CMAuditFormat := recformat< discR, discRp, DegreeFormContent, OddPrimesSubsumed >;

// Record schema for the CM square order construction record.
//   Field, FieldDiscriminant : L and s = D_L = tau^2 (< 0)
//   tau, delta               : trace-zero generator of L, covol(J)/covol(I)
//   Algebra                  : B = (s, -delta)_Q
//   AlgebraDiscriminant      : product of finite ramified primes of B
//   RamifiedPrimes           : finite ramified primes of B
//   Order                    : O (also the first return value)
//   OrderDiscriminant        : Discriminant(O), the reduced discriminant discrd(O)
//   OrderBasis               : the exact quaternion Z-basis [1, r, m1*j, m2*j]
//   PolarGram                : Gram matrix Trd(e_i * conj(e_j)) in OrderBasis
//   Ibasis, Jbasis           : canonical Z-bases of the input lattices I, J
//   Rbasis, Rpbasis, R0basis, Mbasis : canonical Z-bases of R, R', R0, M
//   discR, discRp, discR0    : discriminants of the quadratic orders R, R', R0
//   IsMaximal .. IsBass      : Magma predicates of the order O
//   AuditReasonSets          : the CMAuditFormat record above
//   Fingerprint              : OrderFingerprint(O), the pinned order fingerprint
CMSquareOrderFormat := recformat<
    Field, FieldDiscriminant, tau, delta,
    Algebra, AlgebraDiscriminant, RamifiedPrimes,
    Order, OrderDiscriminant, OrderBasis, PolarGram,
    Ibasis, Jbasis, Rbasis, Rpbasis, R0basis, Mbasis,
    discR, discRp, discR0,
    IsMaximal, IsEichler, IsGorenstein, IsBass,
    AuditReasonSets, Fingerprint >;

// ---------------------------------------------------------------------------
// Internal helpers (file-local).
// ---------------------------------------------------------------------------

// Discriminant of the quadratic order with Z-basis bas (two field elements),
// as the determinant of the field trace form: det(Tr_{L/Q}(b_i b_j)).
function QuadOrderDisc(bas)
    return Determinant(Matrix(Rationals(),
        [[Trace(bas[a]*bas[b]) : b in [1..2]] : a in [1..2]]));
end function;

// Given the canonical Z-basis R0el of an order R0 in L (so 1 in R0), return r
// with R0 = Z*1 + Z*r, r canonically reduced: theta-coordinate > 0 and
// [1,0]-coordinate in [0,1).
function OrderGeneratorR(R0el, L)
    one := L ! 1;
    B := Matrix(Rationals(), [Eltseq(R0el[1]), Eltseq(R0el[2])]);
    v := Vector(Rationals(), Eltseq(one)) * B^-1;   // 1 = v[1]*u1 + v[2]*u2
    assert &and[Denominator(c) eq 1 : c in Eltseq(v)];
    v1 := Integers() ! v[1];
    v2 := Integers() ! v[2];
    g, a, b := Xgcd(v1, v2);
    assert g eq 1;   // 1 is primitive in R0 (an order meets Q in Z)
    // w = [-b, a] completes [v1, v2] to a unimodular basis: v1*a - v2*(-b) = 1.
    r := (-b)*R0el[1] + a*R0el[2];
    rc := Eltseq(r);
    if rc[2] lt 0 then r := -r; rc := Eltseq(r); end if;
    r := r - Floor(rc[1])*one;
    assert CanonicalLatticeBasis([one, r]) eq R0el;   // [1, r] really spans R0
    return r;
end function;

// Do all odd primes dividing n divide D?
function OddPrimesDivide(n, D)
    n := AbsoluteValue(Integers() ! n);
    if n eq 0 then return true; end if;
    D := Integers() ! D;
    return &and[ (IsEven(pe[1]) or D mod pe[1] eq 0) : pe in Factorization(n) ];
end function;

// Core construction. All input validation (same imaginary quadratic field,
// full rank) and the lattice-layer asserts happen inside CMLatticeConstruction,
// which returns field elements in L; L is recovered from that output.
function BuildCMSquareOrder(Ibas, Jbas)
    Rel, Rpel, R0el, Mel, delta := CMLatticeConstruction(Ibas, Jbas);
    L := Universe(Rel);
    Ican := CanonicalLatticeBasis(Ibas);
    Jcan := CanonicalLatticeBasis(Jbas);

    // Deterministic trace-zero generator tau of L with tau^2 = s = D_L < 0.
    OL := MaximalOrder(L);
    s := Discriminant(OL);
    ok, tau := IsSquare(L ! s);
    assert ok;
    if Eltseq(tau)[2] lt 0 then tau := -tau; end if;   // canonical sign
    assert Trace(tau) eq 0 and tau^2 eq (L ! s);

    // Embedding L -> B, x + y*tau |-> x + y*i, expressed on the power basis:
    // theta = L.1 = p + q*tau, so x = a + b*theta |-> (a + b*p) + (b*q)*i.
    te := Eltseq(tau);
    q := 1 / te[2];
    p := -te[1] / te[2];

    B<iB, jB, kB> := QuaternionAlgebra< Rationals() | s, -delta >;
    EmbedL := func< x | (Eltseq(x)[1] + Eltseq(x)[2]*p)*One(B) + (Eltseq(x)[2]*q)*iB >;

    // Pin Magma's j*i = -i*j convention: j*a = conj(a)*j on the embedded L-basis.
    assert jB*iB eq Conjugate(iB)*jB;
    etheta := EmbedL(L.1);
    assert jB*etheta eq Conjugate(etheta)*jB;

    one := L ! 1;
    r := OrderGeneratorR(R0el, L);
    m1 := Mel[1];
    m2 := Mel[2];

    qbasis := [ One(B), EmbedL(r), EmbedL(m1)*jB, EmbedL(m2)*jB ];
    O := QuaternionOrder(Integers(), qbasis);
    discrd := Discriminant(O);

    // ----- Layer 1 asserts -----
    // (1) order closure: every basis product lies in the Z-span of the basis.
    coO := Matrix(Rationals(), [Eltseq(e) : e in qbasis]);
    coOinv := coO^-1;
    for a in [1..4] do
        for b in [1..4] do
            z := Vector(Rationals(), Eltseq(qbasis[a]*qbasis[b])) * coOinv;
            assert &and[Denominator(c) eq 1 : c in Eltseq(z)];
        end for;
    end for;

    // (2) nrd block formula on the deterministic grid.
    for x0, x1, y0, y1 in [-2..2] do
        av := x0*one + x1*r;
        bv := y0*m1 + y1*m2;
        assert Norm(EmbedL(av) + EmbedL(bv)*jB) eq Norm(av) + delta*Norm(bv);
    end for;

    // (3) det(polar Gram) = discrd^2.
    gram := Matrix(Rationals(), 4, 4,
        [[Trace(qbasis[a]*Conjugate(qbasis[b])) : b in [1..4]] : a in [1..4]]);
    assert Determinant(gram) eq discrd^2;

    // (4) reduced trace and norm of every basis element are rational integers.
    assert &and[ Denominator(Trace(e)) eq 1 and Denominator(Norm(e)) eq 1 : e in qbasis ];

    // ----- Discriminants of R, R', R0 -----
    discR  := QuadOrderDisc(Rel);
    discRp := QuadOrderDisc(Rpel);
    discR0 := QuadOrderDisc(R0el);

    // ----- Audit reason sets (degree-form content via the cross-term form) -----
    qm1 := delta*Norm(m1);
    qm2 := delta*Norm(m2);
    qcr := delta*Norm(m1 + m2) - qm1 - qm2;
    assert &and[Denominator(x) eq 1 : x in [qm1, qm2, qcr]];
    content := GCD([Integers() | qm1, qm2, qcr]);
    subsumed := OddPrimesDivide(discR, discrd)
            and OddPrimesDivide(discRp, discrd)
            and OddPrimesDivide(content, discrd);
    audit := rec< CMAuditFormat |
        discR := discR, discRp := discRp,
        DegreeFormContent := content, OddPrimesSubsumed := subsumed >;

    rec_ := rec< CMSquareOrderFormat |
        Field := L, FieldDiscriminant := s, tau := tau, delta := delta,
        Algebra := B, AlgebraDiscriminant := Discriminant(B),
        RamifiedPrimes := RamifiedPrimes(B),
        Order := O, OrderDiscriminant := discrd, OrderBasis := qbasis, PolarGram := gram,
        Ibasis := Ican, Jbasis := Jcan,
        Rbasis := Rel, Rpbasis := Rpel, R0basis := R0el, Mbasis := Mel,
        discR := discR, discRp := discRp, discR0 := discR0,
        IsMaximal := IsMaximal(O), IsEichler := IsEichler(O),
        IsGorenstein := IsGorenstein(O), IsBass := IsBass(O),
        AuditReasonSets := audit, Fingerprint := OrderFingerprint(O) >;

    return O, rec_;
end function;

// ---------------------------------------------------------------------------
// Exported intrinsics.
// ---------------------------------------------------------------------------

intrinsic OrderFingerprint(O::AlgQuatOrd) -> SeqEnum
{ Pinned fingerprint of the definite quaternion order O: the entries of its
  reduced Gram matrix. Deterministic and invariant under rescaling the input
  lattices that produced O (rescaling gives an isometric order). Shared by the
  wrapper and the table generator so a mismatch cannot silently serve a wrong row. }
    return Eltseq(ReducedGramMatrix(O));
end intrinsic;

intrinsic CMSquareOrder(Ibas::SeqEnum, Jbas::SeqEnum) -> AlgQuatOrd, Rec
{ Build the CM square order O = (R meet R') + M*j in B = (s, -delta)_Q from
  Z-bases of full-rank lattices I, J in a common imaginary quadratic field L.
  Returns O and a construction record (fields documented in order.m): R, R', R0,
  M, delta, s, tau, the predicates IsMaximal/IsEichler/IsGorenstein/IsBass,
  discriminants, algebra ramification, the audit reason sets and the fingerprint. }
    return BuildCMSquareOrder(Ibas, Jbas);
end intrinsic;

intrinsic CMSquareOrder(I::RngOrdIdl, J::RngOrdIdl) -> AlgQuatOrd, Rec
{ CMSquareOrder for I, J given as (fractional) ideals of a possibly non-maximal
  quadratic order: the Z-bases are read off with Basis and the general routine runs. }
    LI := NumberField(Order(I));
    LJ := NumberField(Order(J));
    require LI eq LJ : "I and J must be ideals in the same field";
    return BuildCMSquareOrder([LI ! b : b in Basis(I)], [LI ! b : b in Basis(J)]);
end intrinsic;

intrinsic CMSquareOrder(I::RngOrdFracIdl, J::RngOrdFracIdl) -> AlgQuatOrd, Rec
{ CMSquareOrder for fractional-ideal inputs (as above). }
    LI := NumberField(Order(I));
    LJ := NumberField(Order(J));
    require LI eq LJ : "I and J must be ideals in the same field";
    return BuildCMSquareOrder([LI ! b : b in Basis(I)], [LI ! b : b in Basis(J)]);
end intrinsic;
