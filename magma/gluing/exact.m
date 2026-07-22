/*
 * Exact Galois-stable gluing count: the completeness certificate.
 *
 * GaloisStableGluings(E1, E2, ell) returns the EXACT number of Galois-stable
 * anti-symplectic graphs in E1[ell] x E2[ell]. This is the ground truth the
 * analytic pipeline (gluings.m) is checked against: the analytic side finds
 * rational quotients numerically and could in principle undercount (a quotient
 * lost to precision, a conjugate not recognized); the exact count catches that.
 *
 * The unit (quotients). The comparison unit is QUOTIENTS: the RngIntElt returned
 * is the number of M/-M orbits of stable graphs (graph count for ell = 2, half of
 * it for odd ell), and the raw graph count is also recorded. Why not raw graphs:
 * the analytic side recognizes a quotient by its rational moduli, and two graphs
 * that are NOT both Galois-stable can still share moduli (an M whose quotient is
 * isomorphic to a stable one's, e.g. when the target has extra automorphisms), so
 * recognition cannot separate stable graphs from such look-alikes and a per-graph
 * analytic count over-reports. Collapsing to distinct quotients (moduli classes)
 * on both sides removes that ambiguity: M and -M give isomorphic quotients (the -1
 * automorphism of E1), so distinct stable quotients are exactly the M/-M orbits,
 * and the analytic side counts distinct recognized-and-emitted jacobian moduli;
 * products never count toward the certificate.
 *
 * Method (the ell-torsion Galois module is read off a p-adic Galois group; no
 * large number field is ever built, which SplittingField at ell >= 5 makes
 * prohibitive):
 *   1. Short models y^2 = x^3 + A_i x + B_i, Q-isomorphic to E_i (a Q-isomorphism
 *      of elliptic curves is Galois-equivariant on torsion, so E_i[ell] and its
 *      count of stable graphs are unchanged).
 *   2. x-coordinate guard. g_i := DivisionPolynomial(E_i, ell) (roots the
 *      x-coordinates of E_i[ell]); shift g_2 by x -> x - k until coprime to g_1;
 *      G0 := GaloisGroup(g_1 g_2). #G0 = [x-coordinate field : Q]. If it exceeds
 *      DegreeBound, decline (return -1): the Galois group is too big to construct
 *      the module cheaply. This is the cheap early gate (a generic pair at
 *      ell = 7 already blows past 2400).
 *   3. Point polynomial. H_i(T) = prod_{P in E_i[ell]\0}(T - (x(P) + lambda y(P)))
 *      = Res_x(psi_ell(x), (T - x)^2 - lambda^2 (x^3 + A x + B)), lambda chosen so
 *      the ell^2 - 1 values are distinct; its splitting field is the full
 *      ell-torsion field of E_i. For ell = 2 the 2-torsion has y = 0, so H_i is
 *      the 2-division cubic and the points are (x, 0). Shift H_2 to be coprime to
 *      H_1, scale H_1 H_2 to a MONIC INTEGRAL polynomial (roots scaled by the
 *      leading coefficient c) so Magma's GaloisGroup returns its roots
 *      unmodified, and G := GaloisGroup of it. #G = [L : Q], L the compositum of
 *      the two ell-torsion fields. Guard #G against 2 * DegreeBound.
 *   4. Reduce the roots to the residue field F_q of the unramified p-adic
 *      splitting ring, divide by c to recover the u-values, and match each
 *      against E_i(F_q)[ell] (reduction is a Galois-equivariant isomorphism on
 *      the prime-to-p torsion, so the F_ell-module structure over F_q is the one
 *      over Q). Pick bases P_i, Q_i; the discrete logs of the ell^2 points give
 *      rho_i(sigma) in GL2(F_ell) for each generator sigma of G off the
 *      permutation of the torsion points.
 *   5. Anti-symplectic normalization. e_i(P_i, Q_i) = zeta_i, both primitive
 *      ell-th roots of unity in F_q; with zeta_2 = zeta_1^a an anti-symplectic
 *      psi has det = -a^{-1} (mod ell). For ell = 2 the mu_2 pairing is vacuous
 *      (every isomorphism is anti-symplectic), the condition is det = 1.
 *   6. Count M in GL2(F_ell) with det(M) the anti-symplectic value and
 *      rho2(sigma) M = M rho1(sigma) for every generator sigma (equivariance for
 *      generators is equivariance for the whole group, hence full Galois
 *      stability of graph(M)); fold {M, -M} to quotients.
 *
 * The permutation-orientation convention (whether sigma sends root i to i^sigma
 * or i^(sigma^-1)) is irrelevant to the count: the wrong choice replaces each
 * rho_i by its inverse-transpose analogue, and rho2(s) M = M rho1(s) for all s is
 * equivalent to M rho1(s) = rho2(s) M for all s. det(M) uses the pairing, which
 * is read from the points directly, not from the permutation.
 */

// Short Weierstrass A, B (y^2 = x^3 + A x + B) of a curve over Q.
function gsgShortAB(E)
    a := aInvariants(WeierstrassModel(E));
    return a[4], a[5];
end function;

// Point-separating polynomial of E_s : y^2 = x^3 + A x + B at prime ell.
// Odd ell: degree ell^2 - 1, roots x(P) + lambda y(P) over nonzero P in E_s[ell].
// ell = 2: the 2-division cubic (2-torsion has y = 0).
function gsgPointPoly(A, B, ell, lambda)
    if ell eq 2 then
        return DivisionPolynomial(EllipticCurve([A, B]), 2);
    end if;
    P2<xx, TT> := PolynomialRing(Rationals(), 2);
    psi := DivisionPolynomial(EllipticCurve([A, B]), ell);
    return UnivariatePolynomial(Resultant(Evaluate(psi, xx),
        (TT - xx)^2 - lambda^2 * (xx^3 + A*xx + B), xx));
end function;

// Monic integral rescale W = c T (c the leading coefficient) of H: the returned
// polynomial is monic and integral with roots c * (roots of H), so GaloisGroup's
// roots are recoverable (a non-monic H makes GaloisGroup rescale internally and
// GaloisRoot then returns roots of a polynomial that is not H).
function gsgMonicIntegral(H)
    Qx := Parent(H);
    c := LeadingCoefficient(H);
    n := Degree(H);
    return Qx ! [Coefficient(H, i) * c^(n - 1 - i) : i in [0 .. n]], c;
end function;

// u-value -> point map for E_s[ell] over the finite field Fq.
function gsgTorsionMap(A, B, Fq, ell, lambda)
    Eq := EllipticCurve([Fq | A, B]);
    m := AssociativeArray();
    for P in DivisionPoints(Eq ! 0, ell) do
        if P ne Eq ! 0 then
            key := (ell eq 2) select P[1] else P[1] + lambda * P[2];
            m[key] := P;
        end if;
    end for;
    return m, Eq;
end function;

// Basis (P, Q) for one curve's torsion among the root indices idxs, and the
// discrete-log table of the ell^2 combinations a P + b Q. Returns the two chosen
// indices, the two basis points and the table.
function gsgBasis(idxs, pts, ell)
    j1 := idxs[1];
    P := pts[j1];
    j2 := 0;
    for j in idxs do
        if WeilPairing(P, pts[j], ell) ne 1 then j2 := j; break; end if;
    end for;
    error if j2 eq 0, "no second basis point (torsion not full rank)";
    Q := pts[j2];
    tab := AssociativeArray();
    for a in [0 .. ell - 1] do
        for b in [0 .. ell - 1] do
            tab[a*P + b*Q] := <a, b>;
        end for;
    end for;
    return j1, j2, P, Q, tab;
end function;

function gsgFormat()
    return recformat< ell : RngIntElt, group_order : RngIntElt,
                      field_degree : RngIntElt, graph_count : RngIntElt,
                      quotient_count : RngIntElt, matrices : SeqEnum,
                      declined : BoolElt >;
end function;

intrinsic GaloisStableGluings(E1::CrvEll, E2::CrvEll, ell::RngIntElt
    : DegreeBound := 2400) -> RngIntElt, Rec
{The exact number of Galois-stable anti-symplectic gluings of E1 and E2 at prime
ell, counted as QUOTIENTS (M/-M orbits of stable graphs, the unit the analytic
rational-quotient count matches), or -1 when the Galois group is too big to
construct (its order exceeds DegreeBound, resp. 2 * DegreeBound for the full
torsion field). The second value is a record with fields ell, group_order (the
x-coordinate Galois group order), field_degree ([L : Q] for the ell-torsion field
L), quotient_count (the returned number), graph_count (raw stable graphs:
quotient_count for ell = 2, twice it for odd ell), matrices (the stable graph
matrices, in the exact layer's torsion bases, which are NOT the analytic period
bases, so only counts are comparable), and declined. Method: read the ell-torsion
Galois module off the p-adic Galois group of a point-separating polynomial; see
the file header.}
    require IsPrime(ell): "ell must be prime";
    require Type(BaseRing(E1)) eq FldRat and Type(BaseRing(E2)) eq FldRat:
        "E1 and E2 must be defined over Q";
    fmt := gsgFormat();
    Qx := PolynomialRing(Rationals());

    A1, B1 := gsgShortAB(E1);
    A2, B2 := gsgShortAB(E2);

    // x-coordinate guard (cheap early decline).
    g1 := DivisionPolynomial(EllipticCurve([A1, B1]), ell);
    g2 := DivisionPolynomial(EllipticCurve([A2, B2]), ell);
    kx := 0; g2s := g2;
    while GCD(g1, g2s) ne 1 do kx +:= 1; g2s := Evaluate(g2, Qx.1 - kx); end while;
    G0 := GaloisGroup(g1 * g2s);
    if #G0 gt DegreeBound then
        vprintf Gluing: "GaloisStableGluings: x-coordinate group order %o > %o, declining\n", #G0, DegreeBound;
        return -1, rec< fmt | ell := ell, group_order := #G0, field_degree := 0,
            graph_count := -1, quotient_count := -1, matrices := [], declined := true >;
    end if;

    // Full ell-torsion field via the point-separating polynomial.
    lambda := 1;
    H1 := gsgPointPoly(A1, B1, ell, lambda);
    H2 := gsgPointPoly(A2, B2, ell, lambda);
    while ell ne 2 and (not IsSquarefree(H1) or not IsSquarefree(H2)) do
        lambda +:= 1;
        H1 := gsgPointPoly(A1, B1, ell, lambda);
        H2 := gsgPointPoly(A2, B2, ell, lambda);
    end while;
    k := 0; H2s := H2;
    while GCD(H1, H2s) ne 1 do k +:= 1; H2s := Evaluate(H2, Qx.1 - k); end while;
    H := H1 * H2s;
    require IsSquarefree(H): "point polynomial not squarefree (choose another lambda)";
    gW, c := gsgMonicIntegral(H);

    G, _, S := GaloisGroup(gW);
    if #G gt 2 * DegreeBound then
        vprintf Gluing: "GaloisStableGluings: torsion field degree %o > %o, declining\n", #G, 2 * DegreeBound;
        return -1, rec< fmt | ell := ell, group_order := #G0, field_degree := #G,
            graph_count := -1, quotient_count := -1, matrices := [], declined := true >;
    end if;
    vprintf Gluing: "GaloisStableGluings: ell=%o, x-group %o, torsion field degree %o\n", ell, #G0, #G;

    R := Parent(GaloisRoot(1, S));
    OR := Integers(R);
    Fq, red := ResidueClassField(OR);
    cinv := (Fq ! c)^-1;
    m1, _ := gsgTorsionMap(A1, B1, Fq, ell, lambda);
    m2, _ := gsgTorsionMap(A2, B2, Fq, ell, lambda);

    idx1 := []; idx2 := [];
    pts := AssociativeArray();
    for i in [1 .. Degree(G)] do
        u := red(OR ! GaloisRoot(i, S)) * cinv;
        if IsDefined(m1, u) then
            pts[i] := m1[u]; Append(~idx1, i);
        elif IsDefined(m2, u - k) then
            pts[i] := m2[u - k]; Append(~idx2, i);
        else
            error "GaloisStableGluings: a root did not match either torsion module (bad reduction prime?)";
        end if;
    end for;
    require #idx1 eq ell^2 - 1 and #idx2 eq ell^2 - 1:
        "torsion recovery incomplete (unexpected residue-field degree)";

    Fl := GF(ell);
    j1a, j1b, P1, Q1, tab1 := gsgBasis(idx1, pts, ell);
    j2a, j2b, P2, Q2, tab2 := gsgBasis(idx2, pts, ell);

    if ell eq 2 then
        targetdet := Fl ! 1;
    else
        z1 := WeilPairing(P1, Q1, ell);
        z2 := WeilPairing(P2, Q2, ell);
        a := 0;
        for e in [1 .. ell - 1] do if z1^e eq z2 then a := e; break; end if; end for;
        error if a eq 0, "Weil pairings not commensurable (torsion basis degenerate)";
        targetdet := -(Fl ! a)^-1;
    end if;

    gens := GeneratorsSequence(G);
    rho1 := []; rho2 := [];
    for g in gens do
        ca := tab1[pts[j1a^g]]; cb := tab1[pts[j1b^g]];
        da := tab2[pts[j2a^g]]; db := tab2[pts[j2b^g]];
        Append(~rho1, Matrix(Fl, 2, 2, [ca[1], cb[1], ca[2], cb[2]]));
        Append(~rho2, Matrix(Fl, 2, 2, [da[1], db[1], da[2], db[2]]));
    end for;

    stable := [];
    for M in GL(2, Fl) do
        Mm := Matrix(M);
        if Determinant(Mm) ne targetdet then continue; end if;
        if forall{i : i in [1 .. #gens] | rho2[i] * Mm eq Mm * rho1[i]} then
            Append(~stable, Mm);
        end if;
    end for;
    orbits := {};
    for M in stable do Include(~orbits, Min([Eltseq(M), Eltseq(-M)])); end for;

    return #orbits, rec< fmt | ell := ell, group_order := #G0, field_degree := #G,
        graph_count := #stable, quotient_count := #orbits, matrices := stable,
        declined := false >;
end intrinsic;

intrinsic GluingCertificateBlock(E1::CrvEll, E2::CrvEll, n::RngIntElt,
    analyticCount::RngIntElt, Strict::BoolElt : DegreeBound := 2400) -> Tup, MonStgElt
{Compare the analytic rational-quotient count for a prime block against the exact
GaloisStableGluings quotient count and return the info block <n, 1, stable,
analytic, certified> together with this block's proof contribution ("count-matched" or
"traces-only"). The label "count-matched" records only that these two counts agree; it
is not a proof that the emitted curve set is exactly the stable set, since the analytic
recognition is heuristic (README.md gives the precise semantics). Both counts are
distinct rational QUOTIENTS (M/-M orbits on the exact side; distinct
recognized-and-emitted jacobian moduli on the analytic side, products never counting
toward the certificate). On a count disagreement raise "gluing certificate mismatch at
ell=<n>: exact <e> vs analytic <a>". When the exact layer declines (Galois group
over the bound): raise an error if Strict (the Proof := true contract), else return
the uncertified block with "traces-only". Callers gate invocation (Proof "Auto"
runs it for prime n <= 13, Proof := true for every prime with Strict, false never).}
    require IsPrime(n): "n must be prime";
    stable, srec := GaloisStableGluings(E1, E2, n : DegreeBound := DegreeBound);
    if stable lt 0 then
        require not Strict:
            Sprintf("gluing certificate unavailable at ell=%o: exact layer declined (Galois group order %o over the degree bound %o)", n, srec`group_order, DegreeBound);
        return <n, 1, -1, analyticCount, false>, "traces-only";
    end if;
    error if stable ne analyticCount,
        Sprintf("gluing certificate mismatch at ell=%o: exact %o vs analytic %o", n, stable, analyticCount);
    return <n, 1, stable, analyticCount, true>, "count-matched";
end intrinsic;
