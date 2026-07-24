// Standalone driver for the isogeny-primes engine test suite
// (magma/IsogenyPrimes.m).
//
// Data-driven: every expectation (Sage oracle values and recorded engine
// outputs) lives in tests/corpus.txt, whose schema is documented in its
// header; tests/_lib.m parses it and supplies the field/curve constructors,
// tests/_corpus.m attaches the engine and loads the corpus, and
// tests/_suite_body.m holds the section procedures and the RunSuite
// dispatcher (also reused eval-safe by the CI entry Tests/isogenyprimes.m).
// Run from the repository root:
//
//     magma -b tests/test_isogenyprimes.m
//
// Optional -b variables:
//     section:=<name>   run one section: golden, gates, branch1, branch2,
//                       cm, congruence, fixtures, regression (default all)
//     cmscope:=0        skip the cm section (the scope the recorded
//                       regression outputs were taken at)
//     spec:=<path>      engine spec to attach; spec:="" is the red-state
//                       check and must exit 1
//
// Success prints one PASS line per selected section and a final
// "ALL SELECTED SECTIONS PASS" (exit 0). A section failure prints
// "SECTION <name>: FAIL" with the offending entry id; an unknown section or a
// vacuous run (no section executed) prints a top-level "SUITE FAILED" line
// carrying no entry id. Either way the driver quits 1.

load "tests/_corpus.m";
load "tests/_suite_body.m";

if not RunSuite(corpus, section, cmscope) then
    quit 1;
end if;
