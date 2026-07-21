SetQuitOnError(true);  // as in verify.m:1; without it magma -b exits 0 even when an assert/error fires
/* Corpus runner. Usage:
 *   magma -b [tags:=tag1,tag2] [line:=<single corpus line>] check_gluing.m
 * Success = prints PASS per entry and "ALL PASS" at the end; any assert failure aborts.
 * Batch: parallel -a data/gluing_corpus.txt magma -b line:={} check_gluing.m
 *
 * oracle:=1 flag: BHLS-only replay mode, calling Genus2Elliptic2/3 +
 * CanonicalGluingList directly instead of Genus2Gluings, and skipping the
 * proof-field assertion (the oracle path has no certificate).
 *
 * data/gluing_corpus.txt's fields are eval'd below (ParseCurve, field, expected)
 * as trusted, repo-controlled input; this runner is not safe to point at an
 * untrusted corpus file.
 */
AttachSpec("spec");
if not assigned tags then tags := ""; end if;
wanted := tags eq "" select false else Split(tags, ",");
// "assigned oracle" only sees a local binding, not the CLI-level global, once
// evaluated inside a procedure body; capture it in a plain global here (reads
// of an outer global, as opposed to "assigned" checks, do cross into
// procedure scope, the same way "wanted" is read inside RunLine below).
oracleMode := assigned oracle;
// algorithm:=Periods (or Algebraic, only meaningful with tags:=bhls2,bhls3 since
// BHLS covers just n in {2, 3}) forces Genus2Gluings' Algorithm parameter; default
// "Auto". Captured as a plain global for the same cross-scope reason as oracleMode
// above.
algMode := assigned algorithm select algorithm else "Auto";

function ParseCurve(K, s)  // K = shared base field (Rationals() or a NumberField); s = eval'able a-invariant list
    a := eval s;
    if Type(K) eq FldRat then return EllipticCurve([Rationals() | x : x in a]); end if;
    return EllipticCurve([K | K!x : x in a]);
end function;

function Order2Subset(expected)  // expected = list of [fcoeffs, hcoeffs] pairs over Q
    Qx := PolynomialRing(Rationals());
    return [pair : pair in expected |
        #GeometricAutomorphismGroup(HyperellipticCurve(Qx!pair[1], Qx!pair[2])) le 2];
end function;

procedure RunLine(L, ~ran, ~skipped)
    // Blank lines, comment lines, and tag-filtered entries never matched the
    // active run, so they touch neither counter. The oracle-mode skip below is
    // the only path that counts as skipped (matched the filters, deliberately
    // not executed); a completed entry does ran +:= 1 after its PASS printf.
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
    if oracleMode and tag notin ["bhls2", "bhls3"] then
        printf "SKIP %o %o %o (oracle mode covers only bhls2/bhls3)\n", tag, e1, e2;
        skipped +:= 1;
        return;
    end if;
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
    // algorithm:=Periods reconstructs only geometric-automorphism-order-2 quotients
    // (no non-quadratic-twist reconstruction, see gluings.m), and an n=2 block is
    // traces-only on that path, so pin the emitted set to the order-2 subset of the
    // expected curves rather than the Auto-mode count/proof contracts.
    if algMode eq "Periods" then
        if expected cmpne [-1] then
            got := [[Coefficients(f), Coefficients(h)] where f, h := HyperellipticPolynomials(c) : c in cs];
            assert got eq Order2Subset(expected);
        end if;
    else
        if count ge 0 then assert #cs eq count; end if;
        if prf ne "-" and not oracleMode then assert info`proof eq prf; end if;
        if expected cmpne [-1] then
            got := [[Coefficients(f), Coefficients(h)] where f, h := HyperellipticPolynomials(c) : c in cs];
            assert got eq expected;
        end if;
    end if;
    printf "PASS %o %o %o n=%o\n", tag, e1, e2, n;
    ran +:= 1;
end procedure;

ran := 0; skipped := 0;
if assigned line then
    RunLine(line, ~ran, ~skipped);
else
    for L in Split(Read("data/gluing_corpus.txt"), "\n") do RunLine(L, ~ran, ~skipped); end for;
end if;
error if ran eq 0, "no corpus entries ran (filter matched nothing)";
printf "ALL PASS (%o ran, %o skipped)\n", ran, skipped;
quit;
