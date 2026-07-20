/* Corpus runner. Usage:
 *   magma -b [tags:=tag1,tag2] [line:=<single corpus line>] check_gluing.m
 * Success = prints PASS per entry and "ALL PASS" at the end; any assert failure aborts.
 * Batch: parallel -a data/gluing_corpus.txt magma -b line:={} check_gluing.m
 *
 * Temporary oracle:=1 flag (until Genus2Gluings exists): calls
 * Genus2Elliptic2/3 + CanonicalGluingList directly instead of Genus2Gluings,
 * and skips the proof-field assertion (the oracle path has no certificate).
 */
AttachSpec("spec");
if not assigned tags then tags := ""; end if;
wanted := tags eq "" select false else Split(tags, ",");
// "assigned oracle" only sees a local binding, not the CLI-level global, once
// evaluated inside a procedure body; capture it in a plain global here (reads
// of an outer global, as opposed to "assigned" checks, do cross into
// procedure scope, the same way "wanted" is read inside RunLine below).
oracleMode := assigned oracle;

function ParseCurve(field, s)  // s = eval'able a-invariant list
    a := eval s;
    if field cmpeq 0 then return EllipticCurve([Rationals() | x : x in a]); end if;
    K := NumberField(PolynomialRing(Rationals())!field);
    return EllipticCurve([K | K!x : x in a]);
end function;

procedure RunLine(L)
    if #L eq 0 or L[1] eq "#" then return; end if;
    parts := Split(L, ":");
    error if #parts ne 8, "bad corpus line", L;
    tag, fld, e1, e2, ns, cnt, prf, crv :=
        Explode([parts[i] : i in [1..8]]);
    if wanted cmpne false and tag notin wanted then return; end if;
    // Oracle mode only knows the n=2/n=3 BHLS formulas (see dispatch below);
    // any other tag either has n outside {2,3} or is a degenerate pair the
    // formulas were never curated for, so skip rather than silently running
    // the wrong algorithm under the requested n's label.
    if oracleMode and tag notin ["bhls2", "bhls3"] then return; end if;
    field := eval fld;
    E1 := ParseCurve(field, e1); E2 := ParseCurve(field, e2);
    n := StringToInteger(ns);
    if oracleMode then
        cs := CanonicalGluingList(n eq 2 select Genus2Elliptic2(E1, E2) else Genus2Elliptic3(E1, E2));
        info := rec<recformat<proof : MonStgElt> | proof := "-">;
    else
        // eval defers resolution of Genus2Gluings to when this line actually
        // runs: Magma resolves every identifier in a procedure body at
        // definition time, even in a dead branch (confirmed empirically), so
        // a direct call here would fail to compile before Task 7 defines it.
        // Same technique as the SiegelReduce shim in shims.m.
        cs, info := eval "Genus2Gluings(E1, E2, n)";
    end if;
    count := StringToInteger(cnt);
    if count ge 0 then assert #cs eq count; end if;
    if prf ne "-" and not oracleMode then assert info`proof eq prf; end if;
    expected := eval crv;
    if expected cmpne [-1] then
        got := [[Coefficients(f), Coefficients(h)] where f, h := HyperellipticPolynomials(c) : c in cs];
        assert got eq expected;
    end if;
    printf "PASS %o %o %o n=%o\n", tag, e1, e2, n;
end procedure;

if assigned line then
    RunLine(line);
else
    for L in Split(Read("data/gluing_corpus.txt"), "\n") do RunLine(L); end for;
end if;
print "ALL PASS";
quit;
