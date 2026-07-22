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

diffLines := [];
for r in corpus do
    if r`ispair then
        K := BuildField(r`field);
        E1 := BuildCurve(K, r`ainvs);
        E2 := BuildCurve(K, r`ainvs2);
        Lc, infoc := CongruencePrimes(E1, E2);
        EmitCong(~diffLines, r`id, Lc, infoc);
    elif (not IsCMEntry(r)) or cmscope ne "0" then
        K := BuildField(r`field);
        E := BuildCurve(K, r`ainvs);
        Li, infoi := IsogenyPrimes(E);
        EmitIso(~diffLines, r`id, Li, infoi);
    end if;
end for;

Sort(~diffLines);
if out ne "" then
    PrintFile(out, Join(diffLines, "\n") : Overwrite := true);
else
    for line in diffLines do printf "%o\n", line; end for;
end if;
printf "DIFFERENTIAL LINES: %o\n", #diffLines;
