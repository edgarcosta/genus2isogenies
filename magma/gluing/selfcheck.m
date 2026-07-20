intrinsic GluingSelfCheck() -> BoolElt
{Package smoke test: shims and basic invariants. Returns true or raises an assertion.}
    CC := ComplexField(60);
    tau := Matrix(CC, 2, 2, [CC.1, 0, 0, 2*CC.1]);  // already reduced, diagonal
    taured, g := SiegelReduce(tau);
    assert Nrows(taured) eq 2;
    ok, q := RecognizeRational(CC!(22/7) + CC!10^-50*CC.1);
    assert ok and q eq 22/7;
    ok, _ := RecognizeRational(CC!Pi(CC));
    assert not ok;

    // antisymplectic.m checks
    sl2 := func<m | m^3 * &*[Rationals() | 1 - 1/p^2 : p in PrimeDivisors(m)]>;
    for m in [2, 3, 4, 5, 6, 8, 9, 12] do
        A := AntiSymplecticIsomorphisms(m : ModMinus := false);
        assert #A eq Integers()!sl2(m);
        assert forall{M : M in A | Determinant(M) eq -1};
    end for;
    A3 := AntiSymplecticIsomorphisms(3);           // ModMinus default
    assert #A3 eq 12;                               // 24 / 2
    M := AntiSymplecticIsomorphisms(2 : ModMinus := false)[1];
    L := LiftAntiSymplectic(M, 2, 1);
    assert #L eq 8 and forall{N : N in L | Determinant(N) eq -1
        and ChangeRing(N, Integers(2)) eq M};
    C := CRTAntiSymplectic([* AntiSymplecticIsomorphisms(2 : ModMinus := false)[1],
                            AntiSymplecticIsomorphisms(3 : ModMinus := false)[1] *], [2, 3]);
    assert BaseRing(C) eq Integers(6) and Determinant(C) eq -1;

    // modulus.m checks
    N, incon := GluingModulus(EllipticCurve("11a1"), EllipticCurve("11a1"));
    assert incon;                                    // identical curves: all differences zero
    N, incon := GluingModulus(EllipticCurve("11a1"), EllipticCurve("14a1"));
    assert not incon and N ge 1;

    // analytic.m checks: 6 anti-symplectic quotients at n = 2 (|SL2(F2)| = 6,
    // no ModMinus collapse at 2), each classified jacobian or product.
    E1 := EllipticCurve("14a1"); E2 := EllipticCurve("46a1");
    qs := GluedPeriodMatrices(E1, E2, 2 : Precision := 80);
    assert #qs eq 6;
    assert forall{q : q in qs | q`type in ["jacobian", "product"]};

    // EllipticPeriodBasis: orientation normalized and conjugation is an involution.
    ws, Mconj := EllipticPeriodBasis(E1, 60);
    assert Im(ws[1] / ws[2]) gt 0;
    assert Mconj^2 eq IdentityMatrix(Integers(), 2);

    // Step 5 (BHLS numeric cross-check): for the first bhls2 corpus pair, the
    // Igusa-Clebsch invariants of every algebraically glued curve must appear,
    // up to weighted-projective normalization, among the jacobian-type
    // analytic quotients (compared to 20 digits).
    B1 := EllipticCurve([0,0,1,0,-7]); B2 := EllipticCurve([0,0,0,0,4]);
    absIC := func< I | [ I[1]^5 / I[4], I[2]^5 / I[4]^2, I[3]^5 / I[4]^3 ] >;
    exactAbs := [absIC(IgusaClebschInvariants(C)) : C in CanonicalGluingList(Genus2Elliptic2(B1, B2))];
    assert #exactAbs gt 0;
    CC := ComplexField(80);
    bqs := GluedPeriodMatrices(B1, B2, 2 : Precision := 80);
    jacAbs := [absIC(q`invariants) : q in bqs | q`type eq "jacobian"];
    assert #jacAbs gt 0;
    for ea in exactAbs do
        assert exists{ja : ja in jacAbs |
            Max([Abs(CC!ea[t] - ja[t]) / (Abs(CC!ea[t]) + 1) : t in [1..3]]) lt 10^(-20)};
    end for;

    return true;
end intrinsic;
