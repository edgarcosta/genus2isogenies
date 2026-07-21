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
// algorithm:=Periods (or Algebraic) forces Genus2Gluings' Algorithm parameter;
// default "Auto". Captured as a plain global for the same cross-scope reason as
// oracleMode above.
algMode := assigned algorithm select algorithm else "Auto";

function ParseCurve(K, s)  // K = shared base field (Rationals() or a NumberField); s = eval'able a-invariant list
    a := eval s;
    if Type(K) eq FldRat then return EllipticCurve([Rationals() | x : x in a]); end if;
    return EllipticCurve([K | K!x : x in a]);
end function;

function HasHighAutomorphism(expected)  // expected = list of [fcoeffs, hcoeffs] pairs over Q
    Qx := PolynomialRing(Rationals());
    return exists{pair : pair in expected |
        #GeometricAutomorphismGroup(HyperellipticCurve(Qx!pair[1], Qx!pair[2])) gt 2};
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
    // Build K once and share it between E1 and E2: two separate NumberField()
    // calls on the same defining polynomial produce fields that are
    // isomorphic but not cmpeq, which fails GluingModulus/Genus2Gluings'
    // "curves must share a base field" require.
    K := field cmpeq 0 select Rationals() else NumberField(PolynomialRing(Rationals())!field);
    E1 := ParseCurve(K, e1); E2 := ParseCurve(K, e2);
    n := StringToInteger(ns);
    // Genus2Gluings supports any n >= 2 over Q (prime powers via lifting, Task 10;
    // composite via CRT of prime-power blocks, Task 11) and, experimentally
    // (Task 13), prime-power levels over a number field (the nf entry).
    expected := eval crv;
    // algorithm:=Periods only reconstructs curves with geometric automorphism group
    // of order 2 (no non-quadratic-twist reconstruction, see gluings.m); skip the
    // rest so the Periods regression stays a meaningful test of the generic case.
    if not oracleMode and algMode eq "Periods" and expected cmpne [-1] and HasHighAutomorphism(expected) then
        printf "SKIP %o %o %o n=%o (bielliptic/high-automorphism, non-quadratic-twist limitation)\n", tag, e1, e2, n;
        return;
    end if;
    if oracleMode then
        cs := CanonicalGluingList(n eq 2 select Genus2Elliptic2(E1, E2) else Genus2Elliptic3(E1, E2));
        info := rec<recformat<proof : MonStgElt> | proof := "-">;
    else
        // eval defers resolution of Genus2Gluings to when this line actually
        // runs: Magma resolves every identifier in a procedure body at
        // definition time, even in a dead branch (confirmed empirically), so
        // a direct call here would fail to compile before Task 7 defines it.
        // Same technique as the SiegelReduce shim in shims.m.
        cs, info := eval "Genus2Gluings(E1, E2, n : Algorithm := algMode)";
    end if;
    count := StringToInteger(cnt);
    if count ge 0 then assert #cs eq count; end if;
    if prf ne "-" and not oracleMode then assert info`proof eq prf; end if;
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
