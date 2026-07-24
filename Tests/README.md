# Tests

Each `.m` file in this directory is a self-contained Magma test: it must
complete without raising an assertion failure or a runtime error (the
runner treats any escaping error as a failure). Discovered by glob and run
individually through `run_tests.m` from the repo root:

    magma -b target:=Tests/<file>.m exitsignal:="" run_tests.m

`config.m` (repo root) is loaded before every test; a future PR adding a
Magma package should `AttachSpec` it there.

`Tests/ci_smoke.m` is infrastructure, not a regression test: it exists so
CI always has at least one file to discover and run.
