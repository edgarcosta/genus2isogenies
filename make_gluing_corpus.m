/* Regenerates candidate data/gluing_corpus.txt entries for the BHLS oracle
 * formulas (tags bhls2, bhls3): scans Cremona curves for pairs with
 * isomorphic 2-torsion (resp. 3-torsion) splitting fields, runs
 * Genus2Elliptic2/Genus2Elliptic3 + CanonicalGluingList, and prints
 * candidate corpus lines to stdout. Curate the desired lines by hand into
 * data/gluing_corpus.txt (format documented in that file's header).
 *
 * magma -b make_gluing_corpus.m
 */
AttachSpec("spec");

function StripWS(s)
    // Sprintf on nested SeqEnum inserts spaces and, past a width threshold,
    // line breaks; corpus lines must be single-line and space-free.
    return &cat[c : c in Eltseq(s) | c notin {" ", "\n", "\t", "\r"}];
end function;

function CurveList(cs)
    return [[Coefficients(f), Coefficients(h)] where f, h := HyperellipticPolynomials(c) : c in cs];
end function;

// Groups curve indices by the degree of their torsion splitting field
// (cheap to compute for all curves), then only runs the O(group^2)
// IsIsomorphic field comparisons and Genus2Elliptic2/3 attempts within a
// group; matching curves overwhelmingly land in the same degree already,
// so this avoids the full O(#Es^2) field-isomorphism and gluing cost of
// the naive nested loop while returning identical candidate lines.
function Search(Es, torsion, n, tag, target, cap)
    flds := [* torsion(E) : E in Es *];  // List: distinct number fields share no SeqEnum universe
    bydeg := AssociativeArray();
    for i in [1..#Es] do
        d := Degree(flds[i]);
        if IsDefined(bydeg, d) then Append(~bydeg[d], i); else bydeg[d] := [i]; end if;
    end for;
    lines := [];
    used := AssociativeArray();  // per-curve hit count, capped so one curve can't dominate the sample
    for d in Keys(bydeg) do
        idx := bydeg[d];
        for a in [1..#idx] do
            for b in [a+1..#idx] do
                i := idx[a]; j := idx[b];
                ui := IsDefined(used, i) select used[i] else 0;
                uj := IsDefined(used, j) select used[j] else 0;
                if ui ge cap or uj ge cap then continue; end if;
                if not IsIsomorphic(flds[i], flds[j]) then continue; end if;
                E1 := Es[i]; E2 := Es[j];
                try
                    cs := CanonicalGluingList(n eq 2 select Genus2Elliptic2(E1, E2) else Genus2Elliptic3(E1, E2));
                    if #cs gt 0 then
                        crv := CurveList(cs);
                        s := Sprintf("%o:0:%o:%o:%o:%o:-:%o", tag, aInvariants(E1), aInvariants(E2), n, #cs, crv);
                        Append(~lines, StripWS(s));
                        used[i] := ui + 1; used[j] := uj + 1;
                    end if;
                catch e;
                end try;
                if #lines ge target then return lines; end if;
            end for;
        end for;
    end for;
    return lines;
end function;

D := CremonaDatabase();

// 2-torsion block: matching splitting field of the Weierstrass cubic is
// exactly the condition Genus2Elliptic2 itself requires.
Es2 := [EllipticCurve(D, N, i, 1) : i in [1..NumberOfIsogenyClasses(D, N)], N in [11..300]];
torsion2 := func<E | SplittingField(HyperellipticPolynomials(WeierstrassModel(E)))>;
lines2 := Search(Es2, torsion2, 2, "bhls2", 200, 2);

// 3-torsion block: matching splitting field of the 3-division polynomial is
// a necessary (not sufficient) proxy for matching 3-torsion Galois modules;
// Genus2Elliptic3/BHLS3 does the real check and throws or returns [] on most
// candidates, which the try/catch and #cs check skip. Restricted to
// conductor <= 1000 and curves whose division field has degree <= 8 (the
// generic degree is 24; curves that could plausibly match another curve's
// small, special Galois image are the ones worth the O(group^2) search).
Es3 := [EllipticCurve(D, N, i, 1) : i in [1..NumberOfIsogenyClasses(D, N)], N in [11..1000]];
torsion3 := func<E | SplittingField(DivisionPolynomial(E, 3))>;
deg3 := [Degree(torsion3(E)) : E in Es3];
Es3 := [Es3[i] : i in [1..#Es3] | deg3[i] le 8];
lines3 := Search(Es3, torsion3, 3, "bhls3", 60, 2);

print "=== bhls2 ===";
for L in lines2 do print L; end for;
print "=== bhls3 ===";
for L in lines3 do print L; end for;
quit;
