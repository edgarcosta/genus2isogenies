// Loadable preamble for the isogeny-primes test drivers
// (tests/test_isogenyprimes.m and tests/run_differential.m): attach the
// engine, guard its presence, and load the committed data file
// tests/corpus.txt into the `corpus` variable. The parser and the
// field/curve constructors live in the eval-safe tests/_lib.m, loaded here;
// keeping the quit-based guards out of _lib.m lets the CI Tests entry reuse
// the same parser under run_tests.m's eval contract (which forbids
// load/quit). Both drivers load this file first, from the repository root.

load "tests/_lib.m";

// spec:="" is the red-state syntactic check: no engine spec attached; the
// engine-presence guard below then quits 1.
if assigned spec then useSpec := spec; else useSpec := "magma/spec"; end if;
if useSpec ne "" then
    AttachSpec(useSpec);
end if;
if not assigned section then section := "all"; end if;
if not assigned cmscope then cmscope := "1"; end if;
// Engine-presence guard: without the engine intrinsics the driver procedures
// fail at BIND time, which try/catch cannot intercept, so quit with an honest
// verdict before any of that noise.
okEngine, _ := IsIntrinsic("IsogenyPrimes");
okEngine2, _ := IsIntrinsic("CongruencePrimes");
if not (okEngine and okEngine2) then
    printf "SUITE FAILED: engine intrinsics not attached (red state)\n";
    quit 1;
end if;

corpusLoaded := false;
try
    corpus := LoadCorpus("tests/corpus.txt");
    corpusLoaded := true;
catch e
    printf "SUITE FAILED: cannot load tests/corpus.txt: %o\n", e`Object;
end try;
if not corpusLoaded then
    quit 1;
end if;
