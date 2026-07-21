/*
 * AllGenus2Gluings: Algorithm 3.4 (main.tex, algorithm:split non square) driver.
 * Genus2Gluings(F, F', n) glues one PAIR of elliptic curves at one level n; this
 * driver runs Algorithm 3.4's geometrically-nonisogenous-factors sweep: every
 * pair of curves isogenous to the two inputs, times every level n dividing the
 * congruence modulus GluingModulus(E1, E2).
 *
 * Why the sweep, not just Genus2Gluings(E1, E2, n): a gluing of the class
 * representatives (E1, E2) need not be a gluing of E1 and E2 themselves. The
 * LMFDB witness for 2028.a.64896.1 is a (5,5)-gluing of 26a3 x 78a1, NOT of
 * 26a3 x 78a4 (data/gluing_corpus.txt's lmfdb rows; Jac of the witness is
 * Q-isogenous to both products, but the gluing isogeny only descends for one
 * member of 78a4's class). AllGenus2Gluings recovers it by sweeping
 * cls(26a3) x cls(78a4) x {n : n | GluingModulus(26a3, 78a4), n >= 2}.
 *
 * Scope: only Algorithm 3.4's geometrically-nonisogenous branch (E1, E2 not
 * isogenous). GluingModulus's "inconclusive" flag is the necessary-condition
 * proxy for isogenous inputs: a genuine trace difference at a single good
 * prime already proves non-isogenous (Faltings-Serre), so "inconclusive" fires
 * only when every scanned prime agreed. The isogenous (square) case needs
 * Algorithm 3.4's other branch (a 2-power loop over E ~ E'; corollary:kani non
 * square's uniqueness fails there, see main.tex example:non uniqueness), out
 * of scope here.
 *
 * Kani uniqueness (corollary:kani non square): for a genus-2 Jacobian J with
 * NON-isogenous factors, the data (unordered pair of factors up to
 * isomorphism, level n) realizing J as a gluing is unique. So if this sweep
 * ever emits the same curve (up to Q-isomorphism) from two different
 * (pair, level) sources, that is a correctness bug (in this driver or in
 * Genus2Gluings), not a mathematical possibility: checked below and raised as
 * a hard error rather than silently deduplicated.
 */

// Record formats. Plain function/recformat bindings are file-local in Magma
// (a package file cannot see another file's top-level function), so driver.m
// defines its own rather than reusing gluings.m's gluingInfoFmt.
function allgInfoFmt()
    return recformat< N : RngIntElt, sources : SeqEnum, proof : MonStgElt,
                      classsizes : Tup >;
end function;

function allgSourceFmt()
    return recformat< F1 : SeqEnum, F2 : SeqEnum, n : RngIntElt,
                      count : RngIntElt, proof : MonStgElt >;
end function;

// Canonical representative index per element of cls: ids[i] is the smallest j
// such that cls[j] is Q-isomorphic to cls[i] (so ids[i] eq i at each class's
// first occurrence). IsogenousCurves already returns one curve per isomorphism
// class, so this is normally the identity map; computed with IsIsomorphic (not
// j-invariant, which twists share) so a caller-supplied Classes list with
// redundant or twisted entries still identifies factors correctly.
function allgClassRepIds(cls)
    ids := [Integers() | 0 : c in cls];
    for i in [1 .. #cls] do
        for j in [1 .. i] do
            if IsIsomorphic(cls[i], cls[j]) then ids[i] := j; break; end if;
        end for;
    end for;
    return ids;
end function;

// Unordered-pair-of-factor-classes-and-level source key. Gluing (F, F') and
// (F', F) are the same gluing (the graph of an anti-symplectic psi: E1[n] ->
// E2[n] equals the graph of psi^-1: E2[n] -> E1[n], the same subgroup of
// E1[n] x E2[n]), so the Kani comparison must not care which side of the pair
// a factor's class id came from. idA/idB already live in disjoint namespaces
// here (cls1, cls2 come from non-isogenous curves, so can never share an
// isomorphism class) but the pair is still sorted for a well-defined,
// order-independent key.
function allgSourceKey(idA, idB, n)
    tags := Sort(["A" cat IntegerToString(idA), "B" cat IntegerToString(idB)]);
    return <tags[1], tags[2], n>;
end function;

intrinsic AllGenus2Gluings(E1::CrvEll, E2::CrvEll
    : Bound := 500, Classes := false, Proof := "Auto", TraceBound := 1000) -> SeqEnum, Rec
{Algorithm 3.4 (main.tex algorithm:split non square): every genus-2 curve over Q
gluing some curve isogenous to E1 with some curve isogenous to E2, at every level
n >= 2 dividing N := GluingModulus(E1, E2 : Bound := Bound). Sweeps cls(E1) x
cls(E2) x the divisors n of N with n ge 2, where cls(*) is the full isogeny class
(IsogenousCurves) unless Classes := [cls1, cls2] supplies it directly (required
over a number field: there is no isogeny-class-over-K engine yet, so pass the
classes in explicitly). Requires E1, E2 not isogenous (GluingModulus's
inconclusive flag; the isogenous/square case needs Algorithm 3.4's 2-power-loop
branch, not this driver). Proof and TraceBound pass through to every
Genus2Gluings(F, F', n) call. Every emitted curve is checked against Kani
uniqueness (corollary:kani non square): the same curve arising from two distinct
(pair, level) sources raises "uniqueness violation (Kani): same curve from two
gluing data". Returns the CanonicalGluingList merge over every pair and level,
and a metadata record with fields N, sources (one <F1 ainvs, F2 ainvs, n, count,
proof> record per (pair, level) call), proof ("certified" iff every source
certified, else "traces-only"), and classsizes <#cls1, #cls2>.}
    // GluingModulus first, ahead of resolving Classes: an isogenous-input rejection
    // must fire even over a number field with no Classes supplied, without paying for
    // IsogenousCurves (which the FldRat-only require just below would reject anyway).
    N, incon := GluingModulus(E1, E2 : Bound := Bound);
    require not incon:
        "curves appear isogenous; the square case needs the 2-power loop, not this driver";

    if Classes cmpeq false then
        require Type(BaseRing(E1)) eq FldRat and Type(BaseRing(E2)) eq FldRat:
            "AllGenus2Gluings needs Classes := [cls1, cls2] over a number field; there is no isogeny-class-over-K engine yet";
        cls1 := IsogenousCurves(E1);
        cls2 := IsogenousCurves(E2);
    else
        require #Classes eq 2: "Classes must be [cls1, cls2], the isogeny classes of E1 and E2";
        cls1 := Classes[1]; cls2 := Classes[2];
    end if;

    fmt := allgInfoFmt();
    if N le 1 then
        // N = 1: gcd of good-prime trace differences is 1 (GluingModulus's
        // Bound-limited proxy for Algorithm 3.4 Steps 1-2), so no n >= 2 can
        // divide it: certified empty at every level, no pair/level call needed.
        return CanonicalGluingList([]), rec< fmt | N := N, sources := [],
            proof := "certified", classsizes := <#cls1, #cls2> >;
    end if;

    srcfmt := allgSourceFmt();
    ids1 := allgClassRepIds(cls1);
    ids2 := allgClassRepIds(cls2);
    ns := [d : d in Divisors(N) | d ge 2];
    vprintf Gluing: "AllGenus2Gluings: N=%o, %o level(s), sweeping %o x %o pair(s), %o total call(s)\n",
        N, #ns, #cls1, #cls2, #cls1 * #cls2 * #ns;

    allCurves := [];
    sources := [];
    registry := [];   // <curve, sourceKey>: the Kani uniqueness witness set
    for i in [1 .. #cls1] do
        F := cls1[i];
        for j in [1 .. #cls2] do
            Fp := cls2[j];
            for n in ns do
                cs, info := Genus2Gluings(F, Fp, n : Proof := Proof, TraceBound := TraceBound);
                Append(~sources, rec< srcfmt | F1 := aInvariants(F), F2 := aInvariants(Fp),
                    n := n, count := #cs, proof := info`proof >);
                key := allgSourceKey(ids1[i], ids2[j], n);
                for c in cs do
                    matched := false;
                    for r in registry do
                        if IsIsomorphic(c, r[1]) then
                            matched := true;
                            error if r[2] ne key,
                                "uniqueness violation (Kani): same curve from two gluing data";
                            break;
                        end if;
                    end for;
                    if not matched then Append(~registry, <c, key>); end if;
                    Append(~allCurves, c);
                end for;
            end for;
        end for;
    end for;

    proof := (forall{s : s in sources | s`proof eq "certified"}) select "certified" else "traces-only";
    return CanonicalGluingList(allCurves), rec< fmt | N := N, sources := sources,
        proof := proof, classsizes := <#cls1, #cls2> >;
end intrinsic;
