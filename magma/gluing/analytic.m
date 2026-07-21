/*
 * Analytic core of elliptic-curve gluing at a prime level n.
 *
 * Given E1, E2 over Q and a prime n, we enumerate the anti-symplectic
 * isomorphisms psi : E1[n] -> E2[n] (matrices of determinant -1 over
 * Integers(n) in symplectic period bases) and, for each, realise the
 * principally polarized quotient (E1 x E2)/graph(psi) as a period matrix.
 *
 * Analytic picture. E = C / (Z w1 + Z w2) with Im(w1/w2) > 0; the n-torsion
 * point a (w1/n) + b (w2/n) has coordinate vector (a, b). The graph of psi is
 * a maximal isotropic order-n^2 subgroup of (E1 x E2)[n]; its preimage lattice
 * L = Z^4 + Z<(v, psi v)/n> carries a principal polarization, namely n times
 * the product form J = diag(S, S), S = [[0,1],[-1,0]], which is integral and
 * unimodular on L (isotropy of the graph, Pfaffian argument).
 * FrobeniusFormAlternating turns that integral form into the standard genus-2
 * symplectic form and hands back the symplectic basis change; from there the
 * 2x4 big period matrix and the 2x2 small tau follow.
 *
 * Conventions verified empirically, not trusted on paper (see selfcheck.m):
 *   - Periods(E) returns [w1, w2] with w1 real; Im(w1/w2) < 0, so we swap.
 *   - FrobeniusFormAlternating(diag(S,S)) is exactly [[0,I],[-I,0]]; we assert
 *     the transformed form equals that standard genus-2 alternating form.
 *   - tau = A^-1 B must have symmetric, positive-definite-imaginary output; if
 *     the A|B block order gives a non-PD imaginary part we swap the roles.
 *   - Theta(char, 0, tau) with char = [a1, a2, b1, b2]^T (first two a, last two
 *     b) matches the genus-1 product at a diagonal tau to full precision.
 */

// Local height used only for the default-precision heuristic: the largest
// absolute numerator/denominator among the curve's a-invariants.
function CurveHeight(E)
    h := 1;
    for a in aInvariants(E) do
        case Type(a):
            when RngIntElt:
                h := Max(h, Abs(a));
            when FldRatElt:
                h := Max(h, Max(Abs(Numerator(a)), Abs(Denominator(a))));
            else
                for c in Eltseq(a) do
                    cc := Rationals()!c;
                    h := Max(h, Max(Abs(Numerator(cc)), Abs(Denominator(cc))));
                end for;
        end case;
    end for;
    return h;
end function;

// The 10 even genus-2 theta characteristics [a1, a2, b1, b2], a_i, b_i in
// {0, 1/2}: even iff 4 (a . b) is even, i.e. the number of indices with
// a_i = b_i = 1/2 is even (0 or 2).
function EvenThetaChars()
    hs := [0, 1/2];
    out := [];
    for a1 in hs do for a2 in hs do for b1 in hs do for b2 in hs do
        if IsEven(Integers()!(4*(a1*b1 + a2*b2))) then
            Append(~out, [a1, a2, b1, b2]);
        end if;
    end for; end for; end for; end for;
    return out;
end function;

// Is the imaginary part of the 2x2 tau positive definite (leading minors)?
function ImIsPD(tau)
    RR := RealField(Floor(Precision(BaseRing(tau))));
    im := Matrix(RR, 2, 2, [Im(x) : x in Eltseq(tau)]);
    return im[1, 1] gt 0 and Determinant(im) gt 0;
end function;

intrinsic EllipticPeriodBasis(E::CrvEll, prec::RngIntElt) -> SeqEnum, Mtrx
{Period basis [w1, w2] of E to prec decimal digits, normalized so Im(w1/w2) > 0,
together with the 2x2 integer matrix of complex conjugation on that basis
(Conjugate(w_i) = sum_j M[i,j] w_j). Asserts M reproduces the conjugates to half
working precision and is an involution.

Base field. Over Q this is Periods(E). Over a number field K (experimental,
Task 13) it is the period basis of E embedded at ONE fixed real place: the first
InfinitePlaces(K) entry, required real (a complex place would need conjugate-pair
handling, deferred). Periods(E : Precision) over K returns one period pair per
infinite place in InfinitePlaces(K) order; we take place 1's pair. The embedded
curve is defined over R, so complex conjugation acts on its lattice exactly as
over Q and the same solve-and-round recovers M. Both curves in a gluing must use
this same embedding: the complex torus is embedding-dependent, but the rational
gluing output it recognizes is not.}
    require prec ge 10: "precision too low";
    K := BaseRing(E);
    if Type(K) eq FldRat then
        ws := Periods(E : Precision := prec);
    else
        require ISA(Type(K), FldNum):
            "EllipticPeriodBasis: base field must be Q or a number field";
        require IsReal(InfinitePlaces(K)[1]):
            "EllipticPeriodBasis over a number field requires the first infinite place to be real (complex embeddings need conjugate-pair handling, deferred)";
        ws := Periods(E : Precision := prec)[1];   // period pair at InfinitePlaces(K)[1]
    end if;
    if Im(ws[1] / ws[2]) lt 0 then ws := [ws[2], ws[1]]; end if;
    w1 := ws[1]; w2 := ws[2];
    RR := RealField(prec);
    // Solve, per row i, the real 2x2 system R [M[i,1]; M[i,2]] = [Re,Im](conj w_i)
    // where R = [[Re w1, Re w2],[Im w1, Im w2]].
    R := Matrix(RR, 2, 2, [Re(w1), Re(w2), Im(w1), Im(w2)]);
    ent := [];
    for wi in ws do
        x := R^-1 * Matrix(RR, 2, 1, [Re(Conjugate(wi)), Im(Conjugate(wi))]);
        Append(~ent, Round(x[1, 1])); Append(~ent, Round(x[2, 1]));
    end for;
    M := Matrix(Integers(), 2, 2, ent);
    tol := 10^(-(prec div 2));
    assert Abs(M[1, 1] * w1 + M[1, 2] * w2 - Conjugate(w1)) lt tol;
    assert Abs(M[2, 1] * w1 + M[2, 2] * w2 - Conjugate(w2)) lt tol;
    assert M^2 eq IdentityMatrix(Integers(), 2);
    return ws, M;
end intrinsic;

intrinsic GluedBigPeriodMatrix(ws1::SeqEnum, ws2::SeqEnum, psi::Mtrx, n::RngIntElt) -> Mtrx
{The 2x4 big period matrix of (E1 x E2)/graph(psi) with respect to a symplectic
basis for the induced principal polarization. ws1, ws2 are the period bases of
E1, E2 (from EllipticPeriodBasis); psi is a determinant -1 matrix over
Integers(n).}
    require Nrows(psi) eq 2 and Ncols(psi) eq 2: "psi must be 2x2";
    CC := Universe(ws1);
    S := Matrix(Integers(), 2, 2, [0, 1, -1, 0]);
    J := DiagonalJoin(S, S);
    // Graph generators (v, psi v)/n for v = e1, e2 (psi v = column v of psi),
    // scaled by n so n*L = < n e_i, (v, psi v) > is an integer lattice.
    ngens := [
        [Integers()| 1, 0, Integers()!psi[1, 1], Integers()!psi[2, 1]],
        [Integers()| 0, 1, Integers()!psi[1, 2], Integers()!psi[2, 2]]
    ];
    rows := [ [Integers()| n * (i eq j select 1 else 0) : j in [1..4]] : i in [1..4] ]
            cat ngens;
    H := HermiteForm(Matrix(Integers(), #rows, 4, &cat rows));
    basisRows := [i : i in [1..Nrows(H)] | not IsZero(H[i])];
    require #basisRows eq 4: "graph lattice is not rank 4";
    T := Matrix(Rationals(), 4, 4, &cat [ [H[i, j] / n : j in [1..4]] : i in basisRows ]);
    Fq := n * T * ChangeRing(J, Rationals()) * Transpose(T);
    // Integrality is exactly the isotropy of the graph for n*J.
    require forall{c : c in Eltseq(Fq) | Denominator(c) eq 1}:
        "graph is not isotropic: n*J not integral on L";
    F := ChangeRing(Fq, Integers());
    F0, S2 := FrobeniusFormAlternating(F);
    Estd := Matrix(Integers(), 4, 4, [0,0,1,0, 0,0,0,1, -1,0,0,0, 0,-1,0,0]);
    assert F0 eq Estd;
    B := ChangeRing(S2, Rationals()) * T;
    P := Matrix(CC, 2, 4, [ws1[1], ws1[2], 0, 0, 0, 0, ws2[1], ws2[2]]);
    return P * ChangeRing(Transpose(B), CC);
end intrinsic;

intrinsic SmallFromBig(P::Mtrx) -> Mtrx
{Small period matrix tau from a 2x4 big period matrix P = (A | B): tau = A^-1 B,
or B^-1 A if that gives a non-positive-definite imaginary part. Asserts symmetry
to half precision and positive definiteness of Im(tau).}
    require Nrows(P) eq 2 and Ncols(P) eq 4: "P must be 2x4";
    A := Submatrix(P, 1, 1, 2, 2);
    B := Submatrix(P, 1, 3, 2, 2);
    tau := A^-1 * B;
    if not ImIsPD(tau) then tau := B^-1 * A; end if;
    tol := 10^(-(Floor(Precision(BaseRing(P))) div 2));
    assert Abs(tau[1, 2] - tau[2, 1]) lt tol;
    assert ImIsPD(tau);
    return tau;
end intrinsic;

intrinsic NumericInvariants(tau::Mtrx) -> MonStgElt, SeqEnum
{Classify the ppas with small period matrix tau. Siegel-reduce; if the smallest
of the 10 even theta constants is below 10^(-prec/2) the quotient is a product
of elliptic curves and we return "product", [j1, j2] (j-invariants of the
reduced diagonal, or unevaluated zeros if the reduced tau is not diagonal);
otherwise "jacobian", the numeric Igusa-Clebsch quadruple of
y^2 = x (x-1)(x-e1)(x-e2)(x-e3) from RosenhainInvariants.}
    require Nrows(tau) eq 2 and Ncols(tau) eq 2: "tau must be 2x2";
    CC := BaseRing(tau);
    pp := Floor(Precision(CC));
    thr := 10^(-(pp div 2));
    taured := SiegelReduce(tau);
    z := Matrix(CC, 2, 1, [0, 0]);
    tc := [Theta(Matrix(CC, 4, 1, [CC!c : c in ch]), z, taured) : ch in EvenThetaChars()];
    if Min([Abs(t) : t in tc]) lt thr then
        if Abs(taured[1, 2]) lt thr then
            return "product", [jInvariant(taured[1, 1]), jInvariant(taured[2, 2])];
        end if;
        return "product", [CC | 0, 0];
    end if;
    e := Setseq(RosenhainInvariants(taured));
    assert #e eq 3;
    R<x> := PolynomialRing(CC);
    C := HyperellipticCurve(x * (x - 1) * (x - e[1]) * (x - e[2]) * (x - e[3]));
    return "jacobian", IgusaClebschInvariants(C);
end intrinsic;

intrinsic GluingPrecisionHeuristic(E1::CrvEll, E2::CrvEll, n::RngIntElt) -> RngIntElt
{Default analytic precision (decimal digits) for gluing E1, E2 at level n:
40 + 10 n + 5 ceil(log10(1 + H)), H the largest a-invariant height of E1, E2.
Shared by GluedPeriodMatrices (its Precision default) and Genus2Gluings' pass 2
(gluings.m), which needs the same number ahead of calling GluedPeriodMatrices.}
    H := Max(CurveHeight(E1), CurveHeight(E2));
    return 40 + 10 * n + 5 * Ceiling(Log(10, 1 + H));
end intrinsic;

intrinsic GluedPeriodMatrices(E1::CrvEll, E2::CrvEll, n::RngIntElt : Precision := false, Filter := true) -> SeqEnum
{All quotients (E1 x E2)/graph(psi) over the anti-symplectic psi : E1[n] -> E2[n],
for prime-power n = ell^e. Each entry is a record <psi, P, taured, type, invariants>:
psi the determinant -1 matrix, P the 2x4 big period matrix, taured the Siegel-reduced
small period matrix, type one of "jacobian" or "product", invariants the numeric
Igusa-Clebsch quadruple or pair of j-invariants. Precision defaults to
GluingPrecisionHeuristic(E1, E2, n). Filter (default true) restricts psi to the
conjugation-equivariance necessary condition below; Filter := false enumerates
every anti-symplectic psi (the "all quotients over C" contract). This is the
DIRECT enumeration (every det -1 matrix mod n); Genus2Gluings at n = ell^e with
e >= 2 instead lifts a level-ell survivor set (gluings.m), which is the same
rational-quotient set found far more cheaply.}
    require BaseRing(E1) cmpeq BaseRing(E2): "curves must share a base field";
    require IsPrimePower(n): "GluedPeriodMatrices currently supports prime-power n only";
    if Precision cmpeq false then
        prec := GluingPrecisionHeuristic(E1, E2, n);
    else
        prec := Precision;
    end if;
    vprintf Gluing: "GluedPeriodMatrices: n=%o, precision=%o\n", n, prec;
    ws1, c1 := EllipticPeriodBasis(E1, prec);
    ws2, c2 := EllipticPeriodBasis(E2, prec);
    psisAll := AntiSymplecticIsomorphisms(n : ModMinus := false);
    // Conjugation-equivariance filter, empirically derived (see the commit
    // introducing this parameter, not the naive paper formula): a point of
    // E[n] is a coordinate vector v against the w1,w2 basis, but conjugation
    // is defined on the BASIS (Conjugate(w_i) = sum_j c[i,j] w_j, per
    // EllipticPeriodBasis), so it acts on coordinate vectors by c^TRANSPOSE,
    // not c. The graph of psi is stable under the diagonal conjugation action
    // on E1[n] x E2[n] (necessary for it to be Galois-stable, hence for a
    // Q-rational quotient) iff psi * c1^T = c2^T * psi.
    //
    // Checked against 5 independent corpus pairs (4 order-2 bhls3 + the
    // productive lmfdb n=5 pair): this transposed, plus-sign relation is the
    // only one of the 4 candidates (+-, transposed or not) satisfied by every
    // one of the 10 known-Q-rational psi found by brute force; in particular
    // the untransposed minus-sign guess (the naive "antiholomorphic sign
    // flip") holds on all 4 bhls3 pairs but fails outright on the lmfdb pair.
    // Also checked on the exceptional high-automorphism bhls2 corpus entry 1
    // (aut group order up to 24): 3 psis there land on byte-identical rational
    // invariants (the target's extra automorphisms make non-equivariant psi
    // coincide with an equivariant one in moduli), and this relation is the
    // one that picks out exactly the truly equivariant psi among the 3.
    // For prime n the c_i are always conjugate to diag(1,-1) (order 2, det -1,
    // forced since Im(w1/w2) > 0 rules out c_i = +-I), which makes the filter
    // provably tight for odd n: it cuts the n(n-1)(n+1)-size candidate pool to
    // exactly n - 1 survivors (solutions of s t = k, k a fixed nonzero
    // constant, in Z/n), not just "roughly n^2". At a prime power n = ell^e the
    // relation is still a necessary condition for Galois-stability (the conjugation
    // action is n-agnostic), so this direct enumeration stays correct; it is just no
    // longer proven tight, which only affects how many candidates are computed.
    if Filter then
        R := Integers(n);
        c1T := Transpose(ChangeRing(c1, R));
        c2T := Transpose(ChangeRing(c2, R));
        psis := [psi : psi in psisAll | ChangeRing(psi, R) * c1T eq c2T * ChangeRing(psi, R)];
    else
        psis := psisAll;
    end if;
    vprintf Gluing: "GluedPeriodMatrices: %o candidate psi (%o total, Filter:=%o)\n", #psis, #psisAll, Filter;
    fmt := recformat< psi : Mtrx, P : Mtrx, taured : Mtrx, type : MonStgElt, invariants : SeqEnum >;
    out := [];
    for psi in psis do
        P := GluedBigPeriodMatrix(ws1, ws2, psi, n);
        tau := SmallFromBig(P);
        typ, inv := NumericInvariants(tau);
        Append(~out, rec< fmt | psi := psi, P := P, taured := SiegelReduce(tau),
                                type := typ, invariants := inv >);
        vprintf Gluing: "  psi %o -> %o\n", Eltseq(psi), typ;
    end for;
    return out;
end intrinsic;
