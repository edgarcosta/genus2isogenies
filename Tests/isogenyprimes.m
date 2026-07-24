// CI entry for the isogeny-primes engine test suite, discovered by glob and
// run under run_tests.m's contract: the runner evaluates this file's text
// concatenated with "return true;" as a single expression, catching any
// escaping error as a failure. That contract forbids `load` (a declaration
// error inside eval) and `quit` (it kills the runner and bypasses the
// per-file aggregation), so this entry cannot use the standalone driver
// tests/test_isogenyprimes.m directly.
//
// Instead it reuses the exact same suite source: a single nested eval
// concatenates the eval-safe parser (tests/_lib.m) and the section body
// (tests/_suite_body.m) into one unit (declarations are allowed inside a
// single concatenated eval), loads the corpus, and returns RunSuite's
// verdict. A false verdict is turned into a raised error, which run_tests.m
// records as a normal failure rather than a process kill.
//
// The engine spec is attached by config.m before this file runs; the
// idempotent re-attach here keeps the entry self-sufficient.
AttachSpec("magma/spec");

suiteOK := eval (Read("tests/_lib.m") cat "\n" cat Read("tests/_suite_body.m")
    cat "\ncorpus := LoadCorpus(\"tests/corpus.txt\");\nreturn RunSuite(corpus, \"all\", \"1\");");

error if not suiteOK, "Tests/isogenyprimes.m: isogeny-primes suite reported section failures";
