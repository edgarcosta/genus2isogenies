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
    // Filter := false is the pre-Task-8 "all quotients over C" contract, so
    // this count must hold regardless of the conjugation filter's default.
    E1 := EllipticCurve("14a1"); E2 := EllipticCurve("46a1");
    qs := GluedPeriodMatrices(E1, E2, 2 : Precision := 80, Filter := false);
    assert #qs eq 6;
    assert forall{q : q in qs | q`type in ["jacobian", "product"]};

    // EllipticPeriodBasis: orientation normalized and conjugation is an involution.
    ws, Mconj := EllipticPeriodBasis(E1, 60);
    assert Im(ws[1] / ws[2]) gt 0;
    assert Mconj^2 eq IdentityMatrix(Integers(), 2);

    // Step 5 (BHLS numeric cross-check): for the first bhls2 corpus pair, the
    // Igusa-Clebsch invariants of every algebraically glued curve must appear,
    // up to weighted-projective normalization, among the jacobian-type
    // analytic quotients (compared to 20 digits). Filter := false: this
    // validates the analytic engine itself against the algebraic ground
    // truth, a check that predates and is orthogonal to the Task 8 filter.
    B1 := EllipticCurve([0,0,1,0,-7]); B2 := EllipticCurve([0,0,0,0,4]);
    absIC := func< I | [ I[1]^5 / I[4], I[2]^5 / I[4]^2, I[3]^5 / I[4]^3 ] >;
    exactAbs := [absIC(IgusaClebschInvariants(C)) : C in CanonicalGluingList(Genus2Elliptic2(B1, B2))];
    assert #exactAbs gt 0;
    CC := ComplexField(80);
    bqs := GluedPeriodMatrices(B1, B2, 2 : Precision := 80, Filter := false);
    jacAbs := [absIC(q`invariants) : q in bqs | q`type eq "jacobian"];
    assert #jacAbs gt 0;
    for ea in exactAbs do
        assert exists{ja : ja in jacAbs |
            Max([Abs(CC!ea[t] - ja[t]) / (Abs(CC!ea[t]) + 1) : t in [1..3]]) lt 10^(-20)};
    end for;

    // Step 6 (Task 8 conjugation filter): the graph of psi must be stable
    // under complex conjugation, psi * c1^T eq c2^T * psi (c_i the
    // EllipticPeriodBasis conjugation matrix reduced mod n; see analytic.m for
    // the derivation), for every psi the default filter keeps, and the filter
    // must actually cut the candidate pool (else the check below is vacuous).
    G1 := EllipticCurve([0,1,0,4,4]); G2 := EllipticCurve([1,0,1,-1,-2]);  // bhls3 order-2 pair
    _, gc1 := EllipticPeriodBasis(G1, 60);
    _, gc2 := EllipticPeriodBasis(G2, 60);
    R3 := Integers(3);
    gc1T := Transpose(ChangeRing(gc1, R3)); gc2T := Transpose(ChangeRing(gc2, R3));
    gAll := GluedPeriodMatrices(G1, G2, 3 : Precision := 100, Filter := false);
    gFilt := GluedPeriodMatrices(G1, G2, 3 : Precision := 100);
    assert #gFilt lt #gAll;
    assert forall{r : r in gFilt |
        ChangeRing(r`psi, R3) * gc1T eq gc2T * ChangeRing(r`psi, R3)};

    // Every Q-rational jacobian-type quotient found without the filter must
    // still be present (by its recognized, exact invariants; not necessarily
    // via the same psi, since a target with extra automorphisms can let a
    // non-equivariant psi coincide with an equivariant one in moduli) among
    // the filtered candidates.
    for r in gAll do
        if r`type ne "jacobian" then continue; end if;
        okr, ICQr := RecognizeIgusaClebsch(r`invariants);
        if not okr then continue; end if;
        assert exists{s : s in gFilt | s`type eq "jacobian" and
            (oks and ICQs eq ICQr where oks, ICQs := RecognizeIgusaClebsch(s`invariants))};
    end for;

    // exact.m: the completeness certificate's exact side (GaloisStableGluings).
    // 2-division fields are tiny, so it must not decline; the count is in the
    // quotient unit (M/-M orbits).
    cnt := GaloisStableGluings(EllipticCurve("14a1"), EllipticCurve("46a1"), 2);
    assert cnt ge 0;
    // y^2 = x^3 + 1 and y^2 = x^3 + 4 glue at 3 to a single quotient (graphs psi,
    // -psi; one orbit).
    assert GaloisStableGluings(EllipticCurve([0,0,0,0,1]), EllipticCurve([0,0,0,0,4]), 3) eq 1;
    // The 5-torsion pair 66c3 x 11a3 (5 | GluingModulus, so the congruence fast
    // path does not fire): the exact layer proves no rational (5,5)-gluing.
    assert GaloisStableGluings(EllipticCurve([1,0,0,-10065,-389499]), EllipticCurve([0,-1,1,0,0]), 5) eq 0;

    // A certificate count disagreement is a hard error (negative test): feed
    // GluingCertificateBlock a doctored analytic count for a pair whose exact
    // quotient count is 1 and confirm the mismatch error fires with its message.
    caught := false;
    try
        _ := GluingCertificateBlock(EllipticCurve([0,0,0,0,1]), EllipticCurve([0,0,0,0,4]), 3, 999, false);
    catch e
        caught := Position(Sprint(e`Object), "gluing certificate mismatch") gt 0;
    end try;
    assert caught;

    return true;
end intrinsic;
