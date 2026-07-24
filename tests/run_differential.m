// Differential driver for the isogeny-primes engine: recompute the engine's
// answer for every entry of tests/corpus.txt and print one line per entry in
// the recorded differential format, so outputs of different engine versions
// are directly comparable (and feed generate_corpus.py --add-regressions):
//
//     single curve:  id:Kind:Exact:primes:cmdisc
//     curve pair:    id:Kind:Exact:Stabilized:primes:CertificationMethod
//
// Lines are sorted. Run from the repository root; out:=<path> writes the
// lines to a file instead of stdout:
//
//     magma -b cmscope:=0 out:=/tmp/P.out tests/run_differential.m
//
// cmscope:=0 skips the entries whose Sage CM oracle flagged CM (the scope
// the recorded regression outputs are taken at); the default includes them.

load "tests/_corpus.m";

if not assigned out then out := ""; end if;

function BoolStr(b)
    return b select "1" else "0";
end function;

function NormPrimes(L)
    S := Sort(SetToSequence(SequenceToSet([ Integers() ! p : p in L ])));
    return Join([ IntegerToString(p) : p in S ], ",");
end function;

procedure EmitIso(~lines, id, L, info)
    if info`IsCM then cmd := IntegerToString(info`CMOrderDiscriminant); else cmd := "0"; end if;
    Append(~lines, id cat ":" cat info`Kind cat ":" cat BoolStr(info`Exact) cat ":" cat NormPrimes(L) cat ":" cat cmd);
end procedure;

procedure EmitCong(~lines, id, L, info)
    Append(~lines, id cat ":" cat info`Kind cat ":" cat BoolStr(info`Exact) cat ":" cat BoolStr(info`Stabilized) cat ":" cat NormPrimes(L) cat ":" cat info`CertificationMethod);
end procedure;

// In-scope entries: every pair, plus the single curves the differential covers
// at this cmscope (all singles by default; only the non-CM singles at
// cmscope:=0, the scope the pins are baked at). Each engine call is wrapped so
// one entry's failure neither aborts the run nor silently drops the id: an
// errored entry emits a distinguishable diagnostic and is excluded from the
// pin output, and the run then exits nonzero (an incomplete .out must never be
// baked). The successful-entry lines are byte-identical to a clean run.
diffLines := [];
emittedIds := {};
erroredIds := {};
expectedIds := {};
for r in corpus do
    if not (r`ispair or (not IsCMEntry(r)) or cmscope ne "0") then
        continue;
    end if;
    Include(~expectedIds, r`id);
    try
        K := BuildField(r`field);
        if r`ispair then
            E1 := BuildCurve(K, r`ainvs);
            E2 := BuildCurve(K, r`ainvs2);
            Lc, infoc := CongruencePrimes(E1, E2);
            EmitCong(~diffLines, r`id, Lc, infoc);
        else
            E := BuildCurve(K, r`ainvs);
            Li, infoi := IsogenyPrimes(E);
            EmitIso(~diffLines, r`id, Li, infoi);
        end if;
        Include(~emittedIds, r`id);
    catch e
        Include(~erroredIds, r`id);
        printf "DIFFERENTIAL ERROR %o: %o\n", r`id, e`Object;
    end try;
end for;

Sort(~diffLines);
if out ne "" then
    PrintFile(out, Join(diffLines, "\n") : Overwrite := true);
else
    for line in diffLines do printf "%o\n", line; end for;
end if;
printf "DIFFERENTIAL LINES: %o\n", #diffLines;

// The emitted (successful) id set must equal the in-scope set: any errored or
// otherwise missing entry makes them differ and fails the run nonzero, so the
// caller never bakes a partial .out.
missingIds := expectedIds diff (emittedIds join erroredIds);
if emittedIds ne expectedIds then
    printf "DIFFERENTIAL FAILED: %o in-scope, %o emitted; errored=%o missing=%o\n",
        #expectedIds, #emittedIds,
        Sort(SetToSequence(erroredIds)), Sort(SetToSequence(missingIds));
    quit 1;
end if;
