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

    // recognize.m: RecognizeIgusaClebsch on the projectively-I2 = 0 locus (the
    // weight-0 ratio charts). The regression curve y^2 = x^5 + x^4 + 2x^3 + 4x^2 + x
    // has Magma Igusa-Clebsch quadruple [0, -67584, -15925248, -2595225600] (I2 = 0).
    // Embed it at 60 digits, scale coordinate w by a nonreal lambda^w and add a
    // sub-gate residual on the zeroed I2 (the analytic pipeline's roundoff), and
    // require recognition of a weighted-projectively equivalent rational point: the
    // convention-invariant ratios t1 = I4 I6/I10 = -10368/25 and t2 = I4^5/I10^2 =
    // -130842624/625 must come back exactly. (Coordinatewise recovery is ill-posed:
    // weight-0 ratios cannot restore the arbitrary scaling representative.) This
    // fails at aa8251e, where RecognizeIgusaClebsch's I2 pivot rejects the locus.
    ic2 := [CC | 0, -67584, -15925248, -2595225600];
    lam := 1 + CC.1;                                     // nonreal weighted scaling
    ic2s := [CC | ic2[1]*lam^2 + CC!10^(-50), ic2[2]*lam^4, ic2[3]*lam^6, ic2[4]*lam^10];
    okI2, icqI2 := RecognizeIgusaClebsch(ic2s);
    assert okI2;
    assert icqI2[2]*icqI2[3]/icqI2[4] eq -10368/25;
    assert icqI2[2]^5/icqI2[4]^2 eq -130842624/625;
    assert IgusaClebschNearRational(ic2s);

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

    // NormalizedIgusaInvariants must be model-independent, including
    // non-integral models (regression: denominator clearing dropped v[i]).
    Qx<x> := PolynomialRing(Rationals());
    assert NormalizedIgusaInvariants(x^5 + x/3 + 1, Qx!0)
        eq NormalizedIgusaInvariants(9*x^5 + 3*x + 9, Qx!0);

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

    // gluings.m Task 10 (prime-power levels via lifting). The pp corpus pair 15a7 x 30a8
    // at n = 4 = 2^2 is certified EMPTY by reduction: its prime level ell = 2 has 0
    // Galois-stable gluings, and a stable (4,4)-graph reduces mod 2 to a stable (2,2)-graph,
    // so no (4,4)-gluing exists. Block <2, 2, 0, 0, true>, proof "certified" (the firmed
    // corpus count and the certified-empty-by-reduction upgrade).
    csp, ip := Genus2Gluings(EllipticCurve([1,1,1,-80,242]), EllipticCurve([1,0,1,-454,-544]), 4);
    assert #csp eq 0 and ip`proof eq "certified" and ip`blocks eq [<2, 2, 0, 0, true>];

    // The lift SWEEP itself (not the empty short-circuit): the isogenous self-gluing
    // 20a2 x 20a2 at n = 4 skips the certificate, so gluingLiftSweep runs (level-2
    // survivors lifted to mod 4, conjugation-filtered, recognized at full precision). Its
    // single rational quotient is bielliptic and dropped, so 0 curves emit under an e = 2
    // traces-only block.
    Esq := EllipticCurve([0,1,0,-1,0]);
    csq, iq := Genus2Gluings(Esq, Esq, 4);
    assert #csq eq 0 and iq`proof eq "traces-only";
    blk := iq`blocks[1];
    assert blk[1] eq 2 and blk[2] eq 2 and not blk[5];

    // gluings.m Task 11 (composite n = prod ell_i^e_i via CRT composition of prime-power
    // blocks). Three cases:
    //
    // (a) Certified-empty by an exact block. The composite corpus pair 14a4 x 34a3 at n = 6
    // has 2-block exact Galois-stable count 0, so no (2,2)-graph, hence no (6,6)-graph: the
    // composite is certified empty on that block alone (blocks[1] = <2, 1, 0, 0, true>, the 3
    // block left undetermined by the short-circuit). This is the firmed corpus count.
    csc, ic := Genus2Gluings(EllipticCurve([1,0,1,-1,0]), EllipticCurve([1,0,0,-103,-411]), 6);
    assert #csc eq 0 and ic`proof eq "certified" and ic`blocks[1] eq <2, 1, 0, 0, true>;

    // (b) Certified-empty by the congruence fast path (kept first). 11a1 x 17a1 has
    // GluingModulus 1 (conclusive), so no ell^e-part of 6 divides it and every block is
    // congruence-certified empty: proof "certified", 0 curves, before any exact/analytic work.
    cse, ie := Genus2Gluings(EllipticCurve("11a1"), EllipticCurve("17a1"), 6);
    assert #cse eq 0 and ie`proof eq "certified"
        and forall{b : b in ie`blocks | b[3] eq 0 and b[5]};

    // (c) Productive composition + graph-level certificate. 99a2 x 99b3 glue at n = 6: the
    // 2-block (exact graph 2) and 3-block (exact graph 2) CRT-combine, and the level-n
    // recognition finds exactly the product 2 * 2 = 4 rational graphs (2 M/-M quotients), both
    // reconstructing to Q-curves. proof "certified" (every block certified AND the analytic
    // graph count equals the product of block graph counts, the composite certificate).
    csp2, ip2 := Genus2Gluings(EllipticCurve([1,-1,1,-17,30]), EllipticCurve([1,-1,1,-1319,-18084]), 6);
    assert #csp2 eq 2 and ip2`proof eq "certified"
        and forall{b : b in ip2`blocks | b[5]};

    // driver.m (Task 12, Algorithm 3.4): AllGenus2Gluings sweeps the full isogeny classes of
    // both curves and every divisor n >= 2 of GluingModulus. The lmfdb corpus pair 26a3 x 78a4
    // (data/gluing_corpus.txt) has GluingModulus 5 (a single divisor n = 5), so this is a
    // 3 x 4 x 1 = 12-call sweep. The LMFDB witness curve for 2028.a.64896.1 is a (5,5)-gluing
    // of 26a3 x 78a1 (a different member of 78a4's class), NOT of 26a3 x 78a4 directly; the
    // driver must recover it by sweeping the classes, alongside the direct pair's own
    // certified (5,5)-gluing (the corpus's recorded Genus2Gluings(26a3, 78a4, 5) output).
    t0 := Cputime();
    csAllg, infoAllg := AllGenus2Gluings(EllipticCurve("26a3"), EllipticCurve("78a4"));
    vprintf Gluing: "GluingSelfCheck: AllGenus2Gluings(26a3, 78a4) took %o s\n", Cputime(t0);
    Qxd := PolynomialRing(Rationals());
    witnessC := CanonicalGluingList([HyperellipticCurve(Qxd![5,9,6,6,4,1,1], Qxd![0,1,1])])[1];
    directC := CanonicalGluingList([HyperellipticCurve(Qxd![0,-4,-57,-159,246,-105,14], Qxd![0,1,0,1])])[1];
    assert infoAllg`N eq 5 and infoAllg`classsizes eq <3, 4> and infoAllg`proof eq "certified";
    assert exists{c : c in csAllg | IsIsomorphic(c, witnessC)};
    assert exists{c : c in csAllg | IsIsomorphic(c, directC)};

    // Task 13 (number-field inputs, experimental). K = Q(sqrt5); Enf1 = LMFDB
    // 2.2.5.1-31.1-a1 and its Galois conjugate Enf2 (the nf corpus pair).
    // EllipticPeriodBasis at the fixed first (real) infinite place returns a
    // positively oriented lattice whose tau reproduces the embedded jInvariant
    // (validates the place-1 period pair, ordering included) with an involutive
    // conjugation matrix; GluingModulus over K is unchanged; and the full n = 8
    // pipeline runs (~0.2 s) to its firmed empty, traces-only output (the lift
    // sweep finds no rational-looking survivor at level 2 over K).
    Qxnf<xnf> := PolynomialRing(Rationals());
    Knf<th> := NumberField(xnf^2 - xnf - 1);
    Enf1 := EllipticCurve([Knf | 1, 1 + th, th, th, 0]);
    Enf2 := EllipticCurve([Knf | 1, 2 - th, 1 - th, 1 - th, 0]);
    wsnf, Mnf := EllipticPeriodBasis(Enf1, 60);
    assert Im(wsnf[1] / wsnf[2]) gt 0;
    assert Mnf^2 eq IdentityMatrix(Integers(), 2);
    assert Abs(jInvariant(wsnf[1] / wsnf[2])
        - Evaluate(jInvariant(Enf1), InfinitePlaces(Knf)[1] : Precision := 60)) lt 10^-30;
    Nnf, inconnf := GluingModulus(Enf1, Enf2);
    assert Nnf eq 8 and not inconnf;
    csnf, infonf := Genus2Gluings(Enf1, Enf2, 8);
    assert #csnf eq 0 and infonf`proof eq "traces-only"
        and infonf`blocks eq [<2, 3, -1, 0, false>];

    return true;
end intrinsic;
