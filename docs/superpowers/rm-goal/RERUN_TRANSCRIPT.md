# Wave 3 independent clean-rerun transcript

Result: **PASS** for the frozen public seam.  The full 13-test suite passed,
all required direct CLI probes returned the specified JSON status and exit
code, and two independent `D=60` runs were byte-identical.  The only failed
command below is a recorded validation-harness attempt that omitted `HOME`;
Sage rejected that environment before loading the prototype.  Retaining the
real home directory fixed the harness without changing code or repository
state.

## Reproducibility envelope

Validation timestamp: `2026-07-22T04:37:33Z`

Repository and input identities:

```text
HEAD=1c8d56c606443ba622bab8e096e78d3bf630b870
e4df0891b8de036872652feb3317c0568b357708bc51ede7bdac5fcf9ea6afa9  docs/superpowers/rm-goal/FROZEN_SLICE.md
4df8fab62e9c404e8b088e5da1c0ff19722c47301abe66a14a83ecba70537502  docs/superpowers/rm-goal/prototype.sage
2fab63c82cf6b7563a082a170fc9187230ee5b483df44540f14145ddd50ebf98  docs/superpowers/rm-goal/test_prototype.py
SageMath version 10.8, Release Date: 2025-12-18
```

Every passing command ran from
`/home/edgarcosta/projects/genus2isogenies` with this deliberately minimal
environment:

```text
env -i HOME=/home/edgarcosta PATH=/usr/local/bin:/usr/bin:/bin PYTHONHASHSEED=0 LANG=C.UTF-8 LC_ALL=C.UTF-8
```

The explicit `HOME` is not an application input.  It is required by Sage's
launcher.  No network service, Magma process, user configuration variable, or
repository-local cache was added to the environment.

## Minimal-environment harness check

The first attempt intentionally omitted `HOME`:

```sh
/usr/bin/time -p env -i PATH=/usr/local/bin:/usr/bin:/bin PYTHONHASHSEED=0 LANG=C.UTF-8 LC_ALL=C.UTF-8 /usr/local/bin/sage -python docs/superpowers/rm-goal/test_prototype.py -v
```

Exact result:

```text
Error: environment variable $HOME is not set.
Error setting environment variables by sourcing '/home/sage/sage-10.8/src/bin/sage-env';
possibly contact sage-devel (see http://groups.google.com/group/sage-devel).
real 0.11
user 0.03
sys 0.09
```

Exit code: `1`.  This occurred in the Sage launcher before test discovery, so
it is an environment-harness failure rather than a prototype failure.  The
remaining transcript retains the existing real `HOME` and changes nothing
else.

## Full frozen test suite

Command (prefixed by the minimal environment above):

```sh
/usr/local/bin/sage -python docs/superpowers/rm-goal/test_prototype.py -v
```

Captured result:

```text
exit code:       0
wall time:       23.171685 seconds
stdout bytes:    0
stdout sha256:   e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
stderr bytes:    2042
stderr sha256:   0838ff4672bc7e008699e12c2932abf8ccafd97e9d86e138b967493752d1ee15
suite-reported:  Ran 13 tests in 22.991s / OK
```

The exact verbose test list on stderr was:

```text
test_D12_emits_both_positive_unit_square_classes (__main__.PrototypeCLITests.test_D12_emits_both_positive_unit_square_classes) ... ok
test_D145_selects_only_the_order_two_class_from_ordinary_C4 (__main__.PrototypeCLITests.test_D145_selects_only_the_order_two_class_from_ordinary_C4) ... ok
test_D205_rejects_the_ineligible_ordinary_two_class (__main__.PrototypeCLITests.test_D205_rejects_the_ineligible_ordinary_two_class) ... ok
test_D229_uses_the_kernel_not_the_surjective_image_of_squaring (__main__.PrototypeCLITests.test_D229_uses_the_kernel_not_the_surjective_image_of_squaring) ... ok
test_D40_keeps_the_nontrivial_class_in_the_kernel_of_Sq (__main__.PrototypeCLITests.test_D40_keeps_the_nontrivial_class_in_the_kernel_of_Sq) ... ok
test_D5_emits_the_literal_identity_arithmetic_request (__main__.PrototypeCLITests.test_D5_emits_the_literal_identity_arithmetic_request) ... ok
test_D60_emits_the_literal_cartesian_HM_representatives (__main__.PrototypeCLITests.test_D60_emits_the_literal_cartesian_HM_representatives) ... ok
test_bounded_bruteforce_counts_match_genus_theory_and_narrow_two_torsion (__main__.PrototypeCLITests.test_bounded_bruteforce_counts_match_genus_theory_and_narrow_two_torsion) ... ok
test_invalid_CLI_arguments_have_an_explicit_JSON_status (__main__.PrototypeCLITests.test_invalid_CLI_arguments_have_an_explicit_JSON_status) ... ok
test_nonmaximal_order_is_a_validated_rejection_without_HM_data (__main__.PrototypeCLITests.test_nonmaximal_order_is_a_validated_rejection_without_HM_data) ... ok
test_order_discriminant_is_not_accepted_as_a_field_discriminant (__main__.PrototypeCLITests.test_order_discriminant_is_not_accepted_as_a_field_discriminant) ... ok
test_success_document_exposes_every_mandatory_run_check (__main__.PrototypeCLITests.test_success_document_exposes_every_mandatory_run_check) ... ok
test_success_stdout_is_byte_stable (__main__.PrototypeCLITests.test_success_stdout_is_byte_stable) ... ok

----------------------------------------------------------------------
Ran 13 tests in 22.991s

OK
```

`unittest` writes its verbose progress to stderr; this is test-runner output,
not stderr from the prototype subprocesses.  Those are checked separately
below.

## Direct public-CLI probes

Each probe captured stdout and stderr separately in memory.  Every stdout was
valid UTF-8, parsed as exactly one JSON document, and ended in exactly one
newline.

### Maximal arithmetic success: `D=60, f=1`

Command (prefixed by the minimal environment):

```sh
/usr/local/bin/sage -python docs/superpowers/rm-goal/prototype.sage enumerate --field-discriminant 60 --order-conductor 1
```

Exact process metadata:

```text
exit code:       0
wall time:       0.787000 seconds
stdout bytes:    4557
stdout sha256:   b18e5f9f27a7bff631c4b58856d04d6e65cb5a239c46bcbd788240dc9dba66ab
stderr bytes:    0
stderr sha256:   e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
stderr repr:     ''
```

Salient exact decoded JSON values:

```json
{
  "schema": "rm-hm-arithmetic-v1",
  "status": "OK_ARITHMETIC_ONLY",
  "field": {"defining_polynomial": [-15, 0, 1], "discriminant": 60},
  "ordinary_class_group": {"invariants": [2], "order": 2},
  "narrow_class_group": {"invariants": [2, 2], "order": 4},
  "hm_structure": {
    "eligible_ordinary_class_count": 2,
    "ideal_class_condition": "ker(Sq:Cl(F)->Cl+(F))",
    "positive_unit_class_count": 2,
    "record_count": 4
  }
}
```

The four decoded `alpha` values were exactly
`[[1,0],[4,1],[2,0],[8,2]]`; all 15 top-level run checks and every record
check were `true`.

### Validated nonmaximal-order rejection: `D=5, f=2`

Command:

```sh
/usr/local/bin/sage -python docs/superpowers/rm-goal/prototype.sage enumerate --field-discriminant 5 --order-conductor 2
```

Exact process metadata:

```text
exit code:       2
wall time:       0.689247 seconds
stdout bytes:    1209
stdout sha256:   bf274541913512473ddc95cc464aa15dd3a1ae3505610409b61375962cbead44
stderr bytes:    0
stderr sha256:   e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
stderr repr:     ''
```

Salient exact decoded JSON values:

```json
{
  "status": "UNSUPPORTED_NONMAXIMAL_ORDER",
  "field": {"defining_polynomial": [-1, -1, 1], "discriminant": 5},
  "order": {
    "basis": [[1, 0], [0, 2]],
    "conductor": 2,
    "discriminant": 20,
    "index_in_maximal_order": 2,
    "is_maximal": false
  }
}
```

Neither `hm_ideal_classes` nor `hm_records` was present.

### Invalid field discriminant: `D=20, f=1`

Command:

```sh
/usr/local/bin/sage -python docs/superpowers/rm-goal/prototype.sage enumerate --field-discriminant 20 --order-conductor 1
```

Exact process metadata:

```text
exit code:       2
wall time:       0.620842 seconds
stdout bytes:    661
stdout sha256:   411d307784b7553a9d5e6abad48d10a3a73e435386ca47e6f0d3aefd16654b56
stderr bytes:    0
stderr sha256:   e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
stderr repr:     ''
```

Salient exact decoded JSON values:

```json
{
  "input": {"field_discriminant": 20, "order_conductor": 1},
  "schema": "rm-hm-arithmetic-v1",
  "status": "INVALID_FIELD_DISCRIMINANT"
}
```

Neither `field` nor `hm_records` was present.

### Malformed CLI value

Command:

```sh
/usr/local/bin/sage -python docs/superpowers/rm-goal/prototype.sage enumerate --field-discriminant not-an-integer
```

Exact process metadata:

```text
exit code:       2
wall time:       0.604568 seconds
stdout bytes:    598
stdout sha256:   c51f34a791e9b25f6dfa8c68ce2fc1e2061a3214c3544d3a651be25633fee73e
stderr bytes:    86
stderr sha256:   db09d29d078750571e5cd695c2fd5071dcc12a4fb54d9b78c9415f7358eec488
```

Exact decoded stdout status and exact stderr:

```json
{"schema":"rm-hm-arithmetic-v1","status":"INVALID_ARGUMENTS"}
```

```text
invalid arguments: argument --field-discriminant: invalid int value: 'not-an-integer'
```

The compact JSON block above shows the two salient decoded members; the exact
598-byte stdout is identified by its SHA-256 and also contains the frozen
status vocabulary, surface scope, and theorem status.  It contains no
`hm_records`.

## Independent byte-stability rerun

Two new processes ran the same `D=60, f=1` command sequentially, each under
the minimal environment, without retaining output from the earlier direct
probe:

```text
run 1: exit=0, wall=0.773627s, stdout=4557 bytes, stdout_sha256=b18e5f9f27a7bff631c4b58856d04d6e65cb5a239c46bcbd788240dc9dba66ab, stderr=0 bytes
run 2: exit=0, wall=0.693832s, stdout=4557 bytes, stdout_sha256=b18e5f9f27a7bff631c4b58856d04d6e65cb5a239c46bcbd788240dc9dba66ab, stderr=0 bytes
stdout byte comparison: True
stderr byte comparison: True
```

Both empty stderr streams had SHA-256
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`.

## Repository status

Before creating this transcript, the exact root status was:

```text
## main...origin/main
?? .ipynb_checkpoints/
?? Untitled.ipynb
?? Untitled1.ipynb
?? __pycache__/
?? data/output_g2database_2e20.txt
?? private/
?? tests/
```

These entries pre-existed this validation.  `docs/superpowers/` is locally
excluded by `.git/info/exclude`, so adding this assigned transcript does not
alter the porcelain output.  No prototype, test, paper-worktree, report, or
other artifact was edited during this pass.

## Post-review conformance rerun

Validation timestamp: `2026-07-22T04:54:24Z`

Result: **PASS**.  This rerun validates the review-amended working-tree files,
not merely the earlier committed snapshot.  The repository HEAD and exact
input hashes at launch were:

```text
HEAD=c29edc4e108b28aa5c2b66c9ed862c9409ff7f91
e4df0891b8de036872652feb3317c0568b357708bc51ede7bdac5fcf9ea6afa9  docs/superpowers/rm-goal/FROZEN_SLICE.md
d5c42c48d8fcb075c1151fde6cacbeaef8d5af9aad70e0e051e2d1ab21a4f43f  docs/superpowers/rm-goal/prototype.sage
8d3f1df28cc7ee48a6bb04022e0404685f17918ac4bb98bfece90f4ce421a687  docs/superpowers/rm-goal/test_prototype.py
```

The environment and working directory were exactly the same as in the first
passing run: the minimal `env -i` envelope documented above, including Sage's
required real `HOME`, from the repository root.

### Complete public suite

Command:

```sh
/usr/local/bin/sage -python docs/superpowers/rm-goal/test_prototype.py -v
```

Exact captured result:

```text
exit code:       0
wall time:       23.468124 seconds
stdout bytes:    0
stdout sha256:   e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
stderr bytes:    2042
stderr sha256:   1bfddd8661c11eb8ebcd81d325b4a185ace6984fbdd060c47cd084e3955ea314
suite-reported:  Ran 13 tests in 23.363s / OK
```

All 13 test names listed in the initial run again ended in `ok`.  The stderr
hash differs only because the exact suite-reported elapsed time is embedded
in `unittest` output.

### Clean direct probes

The same four public commands produced these exact process results:

```text
D=60,f=1: exit=0, status=OK_ARITHMETIC_ONLY, wall=0.705690s, stdout=4921 bytes, stdout_sha256=80d0bb87f5c7e6851f4a262ae5de0fb6377a27e7c3f6472f8999fcdc514ab9ba, stderr=0 bytes
D=5,f=2:  exit=2, status=UNSUPPORTED_NONMAXIMAL_ORDER, wall=0.637012s, stdout=995 bytes, stdout_sha256=7914b82e8623a14fc243c1aef2a4d2d0fea2692c51bd2c2f4e5f5f5b8e53c409, stderr=0 bytes
D=20,f=1: exit=2, status=INVALID_FIELD_DISCRIMINANT, wall=0.615523s, stdout=661 bytes, stdout_sha256=411d307784b7553a9d5e6abad48d10a3a73e435386ca47e6f0d3aefd16654b56, stderr=0 bytes
malformed: exit=2, status=INVALID_ARGUMENTS, wall=0.639103s, stdout=598 bytes, stdout_sha256=c51f34a791e9b25f6dfa8c68ce2fc1e2061a3214c3544d3a651be25633fee73e, stderr=86 bytes
```

Every stdout parsed as one JSON document.  The three empty stderr streams had
SHA-256
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`.
The malformed case retained the exact expected diagnostic and stderr hash:

```text
invalid arguments: argument --field-discriminant: invalid int value: 'not-an-integer'
db09d29d078750571e5cd695c2fd5071dcc12a4fb54d9b78c9415f7358eec488
```

The amended `D=60` document still has alphas
`[[1,0],[4,1],[2,0],[8,2]]`, four records, and all run/record checks true.
Each record now embeds its referenced ideal exactly: the first two embed
`ideal-0` with norm `1` and HNF `[[1,0],[0,1]]`; the last two embed `ideal-1`
with norm `2` and HNF `[[1,1],[0,2]]`.  The nonmaximal result contains no HM
records and its theorem boundary is now exactly:

```json
{
  "arithmetic_result": "EXACT_ORDER_REJECTION_ONLY",
  "class_unit_enumeration": "NOT_RUN",
  "geometric_interpretation": "NOT_EVALUATED",
  "quotient_ppas": "NOT_COMPUTED",
  "target_deduplication": "NOT_COMPUTED"
}
```

### Fresh byte-stability pair

Two additional `D=60,f=1` processes, independent of the direct probe above,
gave:

```text
run 1: exit=0, wall=0.695811s, stdout=4921 bytes, stdout_sha256=80d0bb87f5c7e6851f4a262ae5de0fb6377a27e7c3f6472f8999fcdc514ab9ba, LF_count=1, final_byte=b'\n', stderr=0 bytes
run 2: exit=0, wall=0.680172s, stdout=4921 bytes, stdout_sha256=80d0bb87f5c7e6851f4a262ae5de0fb6377a27e7c3f6472f8999fcdc514ab9ba, LF_count=1, final_byte=b'\n', stderr=0 bytes
stdout byte comparison: True
stderr byte comparison: True
```

No implementation, test, paper, report, or validation file other than this
assigned transcript was edited by the conformance rerun.
