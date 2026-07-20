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

    return true;
end intrinsic;
