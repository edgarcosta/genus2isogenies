/*
 * Anti-symplectic isomorphisms psi : E1[n] -> E2[n] in coordinates.
 * With both torsion groups written in symplectic bases for the Weil pairing,
 * psi is anti-symplectic iff its matrix has determinant -1 mod n.
 * The graph of psi is a maximal isotropic subgroup of order n^2, and the
 * quotient of E1 x E2 by it is a principally polarized abelian surface.
 */

intrinsic AntiSymplecticIsomorphisms(n::RngIntElt : ModMinus := true) -> SeqEnum
{All 2x2 matrices over Integers(n) of determinant -1. With ModMinus, one
representative of each pair (M, -M) (the quotients they define are related by
the automorphism -1 of E1 and coincide).}
    require n ge 2: "n must be at least 2";
    R := Integers(n);
    G := SL(2, R);
    M0 := Matrix(R, 2, 2, [1, 0, 0, -1]);
    all := [Matrix(g) * M0 : g in G];
    if not ModMinus or n eq 2 then return all; end if;
    seen := {};
    out := [];
    for M in all do
        key := Min([Eltseq(M), Eltseq(-M)]);
        if key notin seen then Include(~seen, key); Append(~out, M); end if;
    end for;
    return out;
end intrinsic;

intrinsic LiftAntiSymplectic(M::Mtrx, ell::RngIntElt, e::RngIntElt) -> SeqEnum
{All lifts of a determinant -1 matrix over Integers(ell^e) to Integers(ell^(e+1))
with determinant -1. There are exactly ell^3.}
    R1 := Integers(ell^(e + 1));
    Mz := Matrix(Integers(), 2, 2, [Integers()!a : a in Eltseq(M)]);
    out := [];
    for a, b, c, d in [0..ell-1] do
        N := ChangeRing(Mz, R1) + ell^e * Matrix(R1, 2, 2, [a, b, c, d]);
        if Determinant(N) eq -1 then Append(~out, N); end if;
    end for;
    return out;
end intrinsic;

intrinsic CRTAntiSymplectic(mats::List, moduli::SeqEnum) -> Mtrx
{Entrywise CRT of matrices over Integers(m_i) for pairwise coprime moduli.
mats is a List (not a SeqEnum) because its entries live over different rings
Integers(m_i), which have no common SeqEnum universe.}
    n := &*moduli;
    ent := [CRT([Integers()!Eltseq(mats[j])[i] : j in [1..#mats]], moduli) : i in [1..4]];
    return Matrix(Integers(n), 2, 2, ent);
end intrinsic;
