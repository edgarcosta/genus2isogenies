// Test suite for the isogeny-primes engine (magma/IsogenyPrimes.m).
//
// Data-driven: every expectation (Sage oracle values and recorded engine
// outputs) lives in tests/corpus.txt, whose schema is documented in its
// header; tests/_corpus.m parses it and supplies the guard preamble and the
// field/curve constructors. Run from the repository root:
//
//     magma -b tests/test_isogenyprimes.m
//
// Optional -b variables:
//     section:=<name>   run one section: golden, gates, branch1, branch2,
//                       cm, congruence, fixtures, regression (default all)
//     cmscope:=0        skip the cm section (the scope the recorded
//                       regression outputs were taken at)
//     spec:=<path>      engine spec to attach; spec:="" is the red-state
//                       check and must exit 1
//
// Success prints one PASS line per selected section and a final
// "ALL SELECTED SECTIONS PASS" (exit 0); any failure prints the failing
// section with the offending entry id and quits 1.

load "tests/_corpus.m";

// Pins the recorded Frobenius data of the two gate entries over Q(sqrt 29)
// (11a1 base-changed; 7 splits, 3 is inert) and the R_q = B_l relation: it
// must hold at the inert prime and must not hold at the split prime unless
// B_l vanishes. The charpolys come from the Sage golden oracle, so they are
// independent of FrobeniusCharpoly.
procedure Test_golden(corpus)
    ngold := 0;
    for r in corpus do
        if r`id notin {"fixture-gate-inert", "fixture-gate-split"} then
            continue;
        end if;
        ngold +:= 1;
        if not assigned r`golden then
            printf "  SKIP %o: golden dropped\n", r`id;
            continue;
        end if;
        K := BuildField(r`field);
        E := BuildCurve(K, r`ainvs);
        sp := r`golden[1];      // [p, Norm(q), charpoly coefficients...]
        ip := r`golden[2];
        qs := PrimeAbove(K, sp[1], 1);
        qi := PrimeAbove(K, ip[1], 1);
        cs := [Integers() ! c : c in Coefficients(FrobeniusCharpoly(E, qs))];
        error if cs ne [sp[j] : j in [3..#sp]],
            Sprintf("%o: split charpoly %o ne recorded %o", r`id, cs,
                    [sp[j] : j in [3..#sp]]);
        ci := [Integers() ! c : c in Coefficients(FrobeniusCharpoly(E, qi))];
        error if ci ne [ip[j] : j in [3..#ip]],
            Sprintf("%o: inert charpoly %o ne recorded %o", r`id, ci,
                    [ip[j] : j in [3..#ip]]);
        if r`id eq "fixture-gate-inert" then
            error if BillereyRq(E, qi) ne BillereyBl(E, ip[1]),
                r`id cat ": inert gate R_q = B_l violated";
        else
            error if not ((BillereyRq(E, qs) ne BillereyBl(E, sp[1]))
                          or (BillereyBl(E, sp[1]) eq 0)),
                r`id cat ": split gate must not hold";
        end if;
    end for;
    error if ngold eq 0, "no gate entries in the corpus";
    printf "SECTION golden: PASS\n";
end procedure;

// Raw engine gates: pure Magma identities of the exposed intrinsics, no
// Sage oracle, so the pinned integers stay literal here rather than in the
// corpus. The model-invariance block, the bracket/star/k=0 identities, and
// the ord([q]) < h_K deviation are each self-contained; xg from the first
// PolynomialRing is reused throughout.
procedure Test_gates()
    // gate: BillereyBl is a model invariant (review repro: 11a1 / Q(sqrt -5), u = 13)
    K := BuildField([5, 0, 1]);
    E := BuildCurve(K, [ [0], [-1], [1], [-10], [-20] ]);
    Es := BuildCurve(K, [ [0], [-169], [2197], [-285610], [-96536180] ]);
    assert IsIsomorphic(E, Es);            // u = 13 rescale is a Weierstrass isomorphism
    q13 := Decomposition(Integers(K), 13)[1][1];
    b1 := BillereyBl(E, 13);
    b2 := BillereyBl(Es, 13);
    assert b1 ne 0;
    assert b1 eq b2;                       // isomorphism invariance
    assert b2 eq BillereyRq(Es, q13);      // R_q = B_l gate on the inert principal q, both models
    // gate: bracket/Adams (PowerCharacteristicPolynomial) and star
    // product (TensorCharacteristicPolynomial) identities, with the
    // k = 0 term ComposePower(.,0) = P(1) pinned through BillereyBl/Rq;
    // integers recorded from the engine on 2026-07-21.
    Zx<xg> := PolynomialRing(Integers());
    K29g<w29g> := BuildField([-29, 0, 1]);
    OK29g := Integers(K29g);
    assert PowerCharacteristicPolynomial((xg-2)*(xg-3), 5) eq (xg - 2^5)*(xg - 3^5);
    // Inert principal q = (3): the R_q k = 0 factor is Evaluate(P, 1).
    Egate29 := EllipticCurve([ K29g | 0, -1, 1, -10, -20 ]);
    q3g := ideal< OK29g | 3 >;
    assert IsPrime(q3g) and Norm(q3g) eq 9;
    assert FrobeniusCharpoly(Egate29, q3g) eq xg^2 + 5*xg + 9;
    P3g := PowerCharacteristicPolynomial(xg^2 + 5*xg + 9, 12);
    assert P3g eq xg^2 - 781282*xg + 282429536481;
    assert Evaluate(P3g, 1) eq 282428755200;          // k = 0 term is P(1)
    assert Evaluate(P3g, 3^12) eq 149653785600;
    assert BillereyRq(Egate29, q3g) eq 42266532377975685120000;
    assert BillereyBl(Egate29, 3) eq 42266532377975685120000;   // inert q=(l): R_q = B_l
    assert BillereyRq(Egate29, q3g) eq Evaluate(P3g, 1)*Evaluate(P3g, 3^12);
    // Split prime 7: the star product of the two prime-above Adams factors.
    Estar29 := EllipticCurve([ K29g | 0, 0, 0, w29g, 1 ]);
    q7sg := [ z[1] : z in Factorization(ideal< OK29g | 7 >) ];
    assert #q7sg eq 2;
    assert { FrobeniusCharpoly(Estar29, q) : q in q7sg } eq { xg^2 - 3*xg + 7, xg^2 + 4*xg + 7 };
    pcpAg := PowerCharacteristicPolynomial(xg^2 - 3*xg + 7, 12);
    pcpBg := PowerCharacteristicPolynomial(xg^2 + 4*xg + 7, 12);
    assert pcpAg eq xg^2 - 136802*xg + 13841287201;
    assert pcpBg eq xg^2 + 153502*xg + 13841287201;
    S7g := TensorCharacteristicPolynomial(pcpAg, pcpBg);   // the star product
    assert S7g eq xg^4 + 20999380604*xg^3 + 202014649792499760006*xg^2
        + 4023087194343502505106185678204*xg
        + 36703368217294125441230211032033660188801;
    assert BillereyBl(Estar29, 7) eq
        8202408623999718705753864179205757894292526928349291149330886396051384418598649856;
    // k = 0 factor Evaluate(.,1) included: dropping it changes B_l.
    assert BillereyBl(Estar29, 7) eq Evaluate(S7g, 1)*Evaluate(S7g, 7^12);
    // gate: pinned h = ord([q]) deviation; K = Q(sqrt -23), h_K = 3.
    // R_q at the inert principal q = (5) uses h = ord([q]) = 1, not h_K;
    // reverting to h_K = 3 changes this integer. Recorded 2026-07-21.
    K23g := BuildField([23, 0, 1]);
    OK23g := Integers(K23g);
    E23g := EllipticCurve([ K23g | 1, 1, 0, 1, 0 ]);
    Cl23g, mCl23g := ClassGroup(K23g);
    assert #Cl23g eq 3;
    q23g := ideal< OK23g | 5 >;
    assert IsPrime(q23g) and Norm(q23g) eq 25 and IsPrincipal(q23g);
    assert Order(q23g @@ mCl23g) eq 1;              // ord([q]) = 1 < h_K = 3
    assert Order(q23g @@ mCl23g) lt #Cl23g;
    q23npg := Factorization(ideal< OK23g | 2 >)[1][1];   // nonprincipal witness
    assert IsPrime(q23npg) and not IsPrincipal(q23npg) and Order(q23npg @@ mCl23g) eq 3;
    assert BillereyRq(E23g, q23npg) eq 75557342874062106394624000000;  // nonprincipal ord([q]) = 3: kills an h := 1 hardcode
    assert FrobeniusCharpoly(E23g, q23g) eq xg^2 + 6*xg + 25;
    P23g := PowerCharacteristicPolynomial(xg^2 + 6*xg + 25, 12);
    assert P23g eq xg^2 - 64250786*xg + 59604644775390625;
    assert Evaluate(P23g, 1) eq 59604644711139840;
    assert Evaluate(P23g, 5^12) eq 103523062500000000;
    assert BillereyRq(E23g, q23g) eq 6170455359721624102560000000000000;
    assert BillereyBl(E23g, 5) eq 6170455359721624102560000000000000;   // inert q=(5): R_q = B_l
    assert BillereyRq(E23g, q23g) eq BillereyBl(E23g, 5);
    printf "SECTION gates: PASS\n";
end procedure;

// Degree-one entries take the exact rational branch: Kind Finite, Exact set,
// and L equals the Sage construction oracle O(E) (complete by Mazur), with
// every oracle ell admitted by MayBeReducible. An entry with a dropped oE
// component only supports containment: the oracle is then a lower bound.
procedure Test_branch1(corpus)
    n := 0;
    for r in corpus do
        if r`ispair or IsCMEntry(r) or not r`deg1 then continue; end if;
        n +:= 1;
        K := BuildField(r`field);
        E := BuildCurve(K, r`ainvs);
        L, info := IsogenyPrimes(E);
        error if info`Source ne "IsogenyPrimes",
            r`id cat ": Source ne IsogenyPrimes";
        error if info`Kind ne "Finite",
            Sprintf("%o: Kind %o ne Finite", r`id, info`Kind);
        error if not info`Exact,
            r`id cat ": degree-one branch must set Exact";
        if not assigned r`oE then
            printf "  SKIP %o: oE dropped\n", r`id;
        elif HasOEDrop(r) then
            // oracle incomplete for this entry: containment only, never equality
            for ell in r`oE do
                error if not MayBeReducible(ell, L, info),
                    Sprintf("%o: oracle ell=%o not admitted (partial oE)", r`id, ell);
            end for;
        else
            error if IntSet(L) ne IntSet(r`oE),
                Sprintf("%o: L=%o ne oracle O(E)=%o", r`id, IntSet(L), IntSet(r`oE));
            for ell in r`oE do
                error if not MayBeReducible(ell, L, info),
                    Sprintf("%o: oracle ell=%o not admitted", r`id, ell);
            end for;
        end if;
    end for;
    error if n eq 0, "no branch1 entries in the corpus";
    printf "SECTION branch1: PASS\n";
end procedure;

// Number-field entries take the Billerey branch: inexact, so the guarantee
// is containment (every construction-certified oracle ell admitted by
// MayBeReducible). Two by-id specials: fixture-localglobal keeps its
// local-global false positive (7 in L although no global 7-isogeny exists),
// and fixture-ex58 must exhaust the B-phase (B_l = 0 for every admissible l)
// and fall through to the R-phase. The Sage soft oracle (reducible_primes)
// is printed for tightness comparison, never asserted.
procedure Test_branch2(corpus)
    n := 0;
    for r in corpus do
        if r`ispair or IsCMEntry(r) or r`deg1 then continue; end if;
        n +:= 1;
        K := BuildField(r`field);
        E := BuildCurve(K, r`ainvs);
        L, info := IsogenyPrimes(E);
        error if info`Source ne "IsogenyPrimes",
            r`id cat ": Source ne IsogenyPrimes";
        error if info`Exact,
            r`id cat ": Billerey branch must not claim Exact";
        if not assigned r`oE then
            printf "  SKIP %o: oE dropped\n", r`id;
        else
            for ell in r`oE do
                error if not MayBeReducible(ell, L, info),
                    Sprintf("%o: oracle ell=%o not admitted (containment)", r`id, ell);
            end for;
        end if;
        if r`id eq "fixture-localglobal" then
            error if 7 notin IntSet(L),
                r`id cat ": local-global false positive 7 must stay in L";
            if assigned r`oE then
                error if 7 in IntSet(r`oE),
                    r`id cat ": oracle found a global 7-isogeny, fixture broken";
            end if;
        end if;
        if r`id eq "fixture-ex58" then
            error if info`BoundsUsed[2] eq 0,
                r`id cat ": R-phase did not run (B_l = 0 for all l)";
        end if;
        if assigned r`soft then
            printf "  branch2 %o soft=%o L=%o\n", r`id, r`soft,
                Sort(SetToSequence(IntSet(L)));
        else
            printf "  branch2 %o soft=UNAVAILABLE L=%o\n", r`id,
                Sort(SetToSequence(IntSet(L)));
        end if;
    end for;
    error if n eq 0, "no branch2 entries in the corpus";
    printf "SECTION branch2: PASS\n";
end procedure;

// CM entries, classified by the Sage CM oracle. The info record must carry
// the oracle's order/fundamental discriminant, conductor, and in-base-field
// flag; Kind is CMFamily iff sqrt(D_F) lies in K (else Finite), Exact iff
// the degree-one branch handled the curve. Finite pins O(E) subset L;
// CMFamily pins the denotation via MayBeReducible. By-id specials:
// fixture-cm-1728-Qi samples the F-family (split and ramified primes are
// construction-reducible), fixture-cm-nonmax pins the family clause at the
// conductor prime p = 2 | f (excluded unless listed), and fixture-h1-cm1728
// re-runs the containment loop on a non-minimal model.
procedure Test_cm(corpus)
    n := 0;
    for r in corpus do
        if r`ispair or not IsCMEntry(r) then continue; end if;
        n +:= 1;
        K := BuildField(r`field);
        E := BuildCurve(K, r`ainvs);
        L, info := IsogenyPrimes(E);
        error if info`Source ne "IsogenyPrimes",
            r`id cat ": Source ne IsogenyPrimes";
        error if not info`IsCM, r`id cat ": engine missed CM";
        error if info`CMOrderDiscriminant ne r`cm[2],
            Sprintf("%o: CMOrderDiscriminant %o ne %o", r`id,
                    info`CMOrderDiscriminant, r`cm[2]);
        error if info`CMFundamentalDiscriminant ne r`cm[3],
            Sprintf("%o: CMFundamentalDiscriminant %o ne %o", r`id,
                    info`CMFundamentalDiscriminant, r`cm[3]);
        error if info`CMConductor ne r`cm[4],
            Sprintf("%o: CMConductor %o ne %o", r`id, info`CMConductor, r`cm[4]);
        inK := r`cm[5] eq 1;
        error if info`CMInBaseField ne inK,
            Sprintf("%o: CMInBaseField %o ne %o", r`id, info`CMInBaseField, inK);
        kind := inK select "CMFamily" else "Finite";
        error if info`Kind ne kind,
            Sprintf("%o: Kind %o ne %o", r`id, info`Kind, kind);
        error if info`Exact ne r`deg1,
            Sprintf("%o: Exact %o ne %o", r`id, info`Exact, r`deg1);
        if not assigned r`oE then
            printf "  SKIP %o: oE dropped\n", r`id;
        elif kind eq "Finite" then
            error if not IntSet(r`oE) subset IntSet(L),
                Sprintf("%o: O(E)=%o not subset L=%o", r`id, IntSet(r`oE), IntSet(L));
        else
            for ell in r`oE do
                error if not MayBeReducible(ell, L, info),
                    Sprintf("%o: oracle ell=%o not admitted (CMFamily denotation)", r`id, ell);
            end for;
        end if;
        if r`id eq "fixture-cm-1728-Qi" then
            if assigned r`oE then
                error if 13 notin IntSet(r`oE),
                    r`id cat ": 13 must be construction-reducible (F-family)";
            end if;
            error if not HasPrimeIsogeny(E, 5),
                r`id cat ": split sample 5 must be reducible (F subset R(E))";
            error if not HasPrimeIsogeny(E, 2),
                r`id cat ": ramified sample 2 must be reducible";
        end if;
        if r`id eq "fixture-cm-nonmax" then
            error if MayBeReducible(2, L, info) ne (2 in IntSet(L)),
                r`id cat ": family clause must reject p | f unless in L";
        end if;
        if r`id eq "fixture-h1-cm1728" and assigned r`oE then
            // non-minimal-model containment: the h1 fixture's whole point
            for ell in r`oE do
                error if not MayBeReducible(ell, L, info),
                    Sprintf("%o: oracle ell=%o not admitted (non-minimal model)", r`id, ell);
            end for;
        end if;
    end for;
    error if n eq 0, "no cm entries in the corpus";
    printf "SECTION cm: PASS\n";
end procedure;

// The fixture-cong-* pairs, dispatched by id/stratum:
//   fixture-cong-deg1-fldnum   degree-one FldNum pair certified AllPrimes by
//                              transport to Q (IsIsogenous);
//   fixture-cong-twist-fldnum  degree-two same-j twist BFS gate: at reduced
//                              bounds (cap below the first separating norm)
//                              the gcd is 0, the certification BFS runs
//                              (BoundsUsed pinned) and must REFUSE the
//                              same-j-but-non-isomorphic node; default
//                              bounds sample the separator and give Finite;
//   fixture-cong-isogenous     isogenous pairs certify AllPrimes, and
//                              KnownIsogenous := true short-circuits with
//                              CertificationMethod Supplied;
//   fixture-cong-twist         Q-level quadratic twist: Finite with the
//                              oracle gcd primes; at the oracle's low cap
//                              the gcd is 0 and degree one returns Undecided
//                              (same-j must never certify);
//   fixture-cong-finite        non-isogenous pairs pin the oracle gcd
//                              primes exactly (fixture-cong-2 is the ell = 2
//                              congruence: both curves have full rational
//                              2-torsion).
procedure Test_congruence(corpus)
    n := 0;
    for r in corpus do
        if not r`ispair or not HasPrefix(r`stratum, "fixture-cong") then
            continue;
        end if;
        n +:= 1;
        K := BuildField(r`field);
        E1 := BuildCurve(K, r`ainvs);
        E2 := BuildCurve(K, r`ainvs2);
        if r`id eq "fixture-cong-deg1-fldnum" then
            L, info := CongruencePrimes(E1, E2);
            error if AbsoluteDegree(BaseRing(E1)) ne 1,
                r`id cat ": fixture must live on a genuine degree-1 FldNum";
            error if info`Kind ne "AllPrimes",
                Sprintf("%o: Kind %o ne AllPrimes (absolute-degree dispatch)", r`id, info`Kind);
            error if not info`Exact, r`id cat ": certified pair must be Exact";
            error if info`CertificationMethod ne "IsIsogenous",
                Sprintf("%o: CertificationMethod %o ne IsIsogenous", r`id,
                        info`CertificationMethod);
        elif r`id eq "fixture-cong-twist-fldnum" then
            // Reduced bounds keep the BFS cheap; the cap sits below the
            // first separating norm 7 so G = 0 and the certification BFS
            // runs (BoundsUsed[2] ne 0 distinguishes it from the degree-1
            // short circuit). Recorded from the engine on 2026-07-21.
            Ltw, infoTw := CongruencePrimes(E1, E2 : KnownIsogenous := false,
                NormBound := 6, MaxNormBound := 6,
                CertificationPrimeBound := 5, CertificationDepth := 1);
            error if jInvariant(E1) ne jInvariant(E2),
                r`id cat ": fixture pair must share j (in-suite witness)";
            error if infoTw`Kind ne "Undecided",
                Sprintf("%o: Kind %o ne Undecided (BFS must reject the same-j node)", r`id, infoTw`Kind);
            error if infoTw`Kind eq "AllPrimes",
                r`id cat ": twist guard: same j must not certify";
            error if infoTw`BoundsUsed ne [6, 5, 1],
                Sprintf("%o: BoundsUsed %o ne [6,5,1]", r`id, infoTw`BoundsUsed);
            error if infoTw`BoundsUsed[2] eq 0,
                r`id cat ": BFS path expected, not the degree-1 short circuit";
            L, info := CongruencePrimes(E1, E2 : KnownIsogenous := false);
            error if info`Kind ne "Finite",
                Sprintf("%o: default bounds: Kind %o ne Finite (norm-7 separator)", r`id, info`Kind);
            error if 2 notin IntSet(L),
                r`id cat ": quadratic twist must be 2-congruent";
            error if IntSet(L) ne {2},
                Sprintf("%o: default bounds: L=%o ne {2}", r`id, IntSet(L));
        elif r`stratum eq "fixture-cong-isogenous" then
            L, info := CongruencePrimes(E1, E2);
            error if info`Source ne "CongruencePrimes",
                r`id cat ": Source ne CongruencePrimes";
            error if info`Kind ne "AllPrimes",
                Sprintf("%o: Kind %o ne AllPrimes (isogenous over Q)", r`id, info`Kind);
            L2, info2 := CongruencePrimes(E1, E2 : KnownIsogenous := true);
            error if info2`Kind ne "AllPrimes",
                Sprintf("%o: KnownIsogenous: Kind %o ne AllPrimes", r`id, info2`Kind);
            error if info2`CertificationMethod ne "Supplied",
                Sprintf("%o: KnownIsogenous: CertificationMethod %o ne Supplied",
                        r`id, info2`CertificationMethod);
        elif r`stratum eq "fixture-cong-twist" then
            L, info := CongruencePrimes(E1, E2);
            error if info`Kind ne "Finite",
                Sprintf("%o: Kind %o ne Finite", r`id, info`Kind);
            if assigned r`congprimes then
                error if IntSet(L) ne IntSet(r`congprimes),
                    Sprintf("%o: L=%o ne oracle gcd primes %o", r`id, IntSet(L),
                            IntSet(r`congprimes));
            end if;
            // At the oracle's low cap the gcd is formally zero (lowkind
            // ZeroAtCap) and degree one short-circuits to Undecided; a
            // same-j pair must never certify AllPrimes from that state.
            if assigned r`lowkind then
                error if r`lowkind ne "ZeroAtCap",
                    Sprintf("%o: lowbound oracle kind %o ne ZeroAtCap (corpus broken)", r`id, r`lowkind);
                L2, info2 := CongruencePrimes(E1, E2 :
                    NormBound := r`lowbound, MaxNormBound := r`lowbound);
                error if info2`Kind ne "Undecided",
                    Sprintf("%o: low cap: Kind %o ne Undecided", r`id, info2`Kind);
                error if info2`Kind eq "AllPrimes",
                    r`id cat ": twist guard: same-j must not certify";
            end if;
        else   // fixture-cong-finite
            L, info := CongruencePrimes(E1, E2);
            error if info`Kind ne "Finite",
                Sprintf("%o: Kind %o ne Finite", r`id, info`Kind);
            if assigned r`congprimes then
                if r`id eq "fixture-cong-2" and 2 in IntSet(r`congprimes) then
                    error if 2 notin IntSet(L),
                        r`id cat ": ell = 2 congruence missed";
                end if;
                error if IntSet(L) ne IntSet(r`congprimes),
                    Sprintf("%o: L=%o ne oracle gcd primes %o", r`id, IntSet(L),
                            IntSet(r`congprimes));
            end if;
        end if;
    end for;
    error if n eq 0, "no fixture-cong entries in the corpus";
    printf "SECTION congruence: PASS\n";
end procedure;

// fixture-deg1-fldnum: dispatch is on ABSOLUTE degree, not Magma type, so a
// genuine degree-one FldNum must take the exact rational branch.
procedure Test_fixtures(corpus)
    found := false;
    for r in corpus do
        if r`id ne "fixture-deg1-fldnum" then continue; end if;
        found := true;
        K := BuildField(r`field);
        E := BuildCurve(K, r`ainvs);
        error if AbsoluteDegree(BaseRing(E)) ne 1,
            r`id cat ": fixture must live on a genuine degree-1 FldNum";
        L, info := IsogenyPrimes(E);
        error if not info`Exact, r`id cat ": branch-1 semantics expected";
        error if info`Kind ne "Finite",
            Sprintf("%o: Kind %o ne Finite", r`id, info`Kind);
        if assigned r`oE then
            if HasOEDrop(r) then
                for ell in r`oE do
                    error if not MayBeReducible(ell, L, info),
                        Sprintf("%o: oracle ell=%o not admitted (partial oE)", r`id, ell);
                end for;
            else
                error if IntSet(L) ne IntSet(r`oE),
                    Sprintf("%o: L=%o ne oracle O(E)=%o", r`id, IntSet(L), IntSet(r`oE));
            end if;
        end if;
    end for;
    error if not found, "fixture-deg1-fldnum missing from the corpus";
    printf "SECTION fixtures: PASS\n";
end procedure;

// Recorded engine outputs (the corpus reg* fields, baked by
// generate_corpus.py --add-regressions from a differential run): a
// single-model regression pin, not an external oracle. Entries run in id
// order, matching the sorted differential .out file the pins came from.
procedure Test_regression(corpus)
    ids := [];
    byid := AssociativeArray();
    for r in corpus do
        if not assigned r`regkind then continue; end if;
        Append(~ids, r`id);
        byid[r`id] := r;
    end for;
    error if #ids eq 0, "no regression pins in the corpus";
    Sort(~ids);
    for id in ids do
        r := byid[id];
        K := BuildField(r`field);
        if r`ispair then
            E1 := BuildCurve(K, r`ainvs);
            E2 := BuildCurve(K, r`ainvs2);
            L, info := CongruencePrimes(E1, E2);
            error if IntSet(L) ne IntSet(r`regprimes),
                Sprintf("%o: L=%o ne recorded %o", id, IntSet(L), IntSet(r`regprimes));
            error if info`Kind ne r`regkind,
                Sprintf("%o: Kind %o ne recorded %o", id, info`Kind, r`regkind);
            error if info`Exact ne r`regexact,
                Sprintf("%o: Exact %o ne recorded %o", id, info`Exact, r`regexact);
            error if info`Stabilized ne r`regstab,
                Sprintf("%o: Stabilized %o ne recorded %o", id, info`Stabilized, r`regstab);
            error if info`CertificationMethod ne r`regcert,
                Sprintf("%o: CertificationMethod %o ne recorded %o", id,
                        info`CertificationMethod, r`regcert);
        else
            E := BuildCurve(K, r`ainvs);
            L, info := IsogenyPrimes(E);
            error if IntSet(L) ne IntSet(r`regprimes),
                Sprintf("%o: L=%o ne recorded %o", id, IntSet(L), IntSet(r`regprimes));
            error if info`Kind ne r`regkind,
                Sprintf("%o: Kind %o ne recorded %o", id, info`Kind, r`regkind);
            error if info`Exact ne r`regexact,
                Sprintf("%o: Exact %o ne recorded %o", id, info`Exact, r`regexact);
            cmd := info`IsCM select info`CMOrderDiscriminant else 0;
            error if cmd ne r`regcmdisc,
                Sprintf("%o: CM discriminant %o ne recorded %o", id, cmd, r`regcmdisc);
        end if;
    end for;
    printf "SECTION regression: PASS\n";
end procedure;

// Magma -b continues past an uncaught top-level runtime error onto the next
// top-level statement, so an unconditional footer would still print PASS
// after a section blew up: wrap each call so a caught error flips ok. A
// procedure that failed to compile (missing intrinsic) is never bound, and
// calling an unbound name raises an error try/catch does NOT intercept, so
// guard each call with `assigned` first.
ok := true;
if section eq "all" or section eq "golden" then
    if assigned Test_golden then
        try
            Test_golden(corpus);
        catch e
            ok := false;
            printf "SECTION golden: FAIL: %o\n", e`Object;
        end try;
    else
        ok := false;
        printf "SECTION golden: FAIL: procedure not declared (spec/engine not attached?)\n";
    end if;
end if;
if section eq "all" or section eq "gates" then
    if assigned Test_gates then
        try
            Test_gates();
        catch e
            ok := false;
            printf "SECTION gates: FAIL: %o\n", e`Object;
        end try;
    else
        ok := false;
        printf "SECTION gates: FAIL: procedure not declared (spec/engine not attached?)\n";
    end if;
end if;
if section eq "all" or section eq "branch1" then
    if assigned Test_branch1 then
        try
            Test_branch1(corpus);
        catch e
            ok := false;
            printf "SECTION branch1: FAIL: %o\n", e`Object;
        end try;
    else
        ok := false;
        printf "SECTION branch1: FAIL: procedure not declared (spec/engine not attached?)\n";
    end if;
end if;
if section eq "all" or section eq "branch2" then
    if assigned Test_branch2 then
        try
            Test_branch2(corpus);
        catch e
            ok := false;
            printf "SECTION branch2: FAIL: %o\n", e`Object;
        end try;
    else
        ok := false;
        printf "SECTION branch2: FAIL: procedure not declared (spec/engine not attached?)\n";
    end if;
end if;
if (section eq "all" or section eq "cm") and cmscope ne "0" then
    if assigned Test_cm then
        try
            Test_cm(corpus);
        catch e
            ok := false;
            printf "SECTION cm: FAIL: %o\n", e`Object;
        end try;
    else
        ok := false;
        printf "SECTION cm: FAIL: procedure not declared (spec/engine not attached?)\n";
    end if;
end if;
if section eq "all" or section eq "congruence" then
    if assigned Test_congruence then
        try
            Test_congruence(corpus);
        catch e
            ok := false;
            printf "SECTION congruence: FAIL: %o\n", e`Object;
        end try;
    else
        ok := false;
        printf "SECTION congruence: FAIL: procedure not declared (spec/engine not attached?)\n";
    end if;
end if;
if section eq "all" or section eq "fixtures" then
    if assigned Test_fixtures then
        try
            Test_fixtures(corpus);
        catch e
            ok := false;
            printf "SECTION fixtures: FAIL: %o\n", e`Object;
        end try;
    else
        ok := false;
        printf "SECTION fixtures: FAIL: procedure not declared (spec/engine not attached?)\n";
    end if;
end if;
// cmscope does NOT gate this section: every regression entry is
// shared-scope by construction (the recorded differential runs at
// cmscope:=0, so no CM-only entry ever gets a pin).
if section eq "all" or section eq "regression" then
    if assigned Test_regression then
        try
            Test_regression(corpus);
        catch e
            ok := false;
            printf "SECTION regression: FAIL: %o\n", e`Object;
        end try;
    else
        ok := false;
        printf "SECTION regression: FAIL: procedure not declared (spec/engine not attached?)\n";
    end if;
end if;
// Failure quits nonzero so batch drivers and CI see it; success falls
// through to PASS and Magma's default exit 0.
if not ok then
    printf "SUITE FAILED\n";
    quit 1;
end if;
printf "ALL SELECTED SECTIONS PASS\n";
