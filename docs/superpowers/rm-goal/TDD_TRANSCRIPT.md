# Maximal-RM prototype TDD transcript

All implementation tests exercise only the frozen process boundary:

```sh
PYTHONHASHSEED=0 sage -python docs/superpowers/rm-goal/test_prototype.py -v
```

No mock or private-function test is used.  The exact RED excerpts below were
captured while the public seam was developed.  A documentation rewrite was
interrupted after deleting the first draft of this file; this reconstruction
uses those already-captured terminal excerpts and the final independent GREEN
run.  No missing intermediate output is invented.

## Tracer cycle 1: the identity request for `D=5`

The first test was run before `prototype.sage` existed.  It failed with exit
1 while trying to decode the empty stdout:

```text
test_D5_emits_the_literal_identity_arithmetic_request (__main__.PrototypeCLITests.test_D5_emits_the_literal_identity_arithmetic_request) ... ERROR

======================================================================
ERROR: test_D5_emits_the_literal_identity_arithmetic_request (__main__.PrototypeCLITests.test_D5_emits_the_literal_identity_arithmetic_request)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/edgarcosta/projects/genus2isogenies/docs/superpowers/rm-goal/test_prototype.py", line 41, in test_D5_emits_the_literal_identity_arithmetic_request
    completed, document = run_enumerate(5)
  File "/home/edgarcosta/projects/genus2isogenies/docs/superpowers/rm-goal/test_prototype.py", line 36, in run_enumerate
    return completed, json.loads(completed.stdout)
json.decoder.JSONDecodeError: Expecting value: line 1 column 1 (char 0)

----------------------------------------------------------------------
Ran 1 test in 0.046s

FAILED (errors=1)
```

After the minimal exact field/order/unit path was added, the same command was
GREEN:

```text
test_D5_emits_the_literal_identity_arithmetic_request (__main__.PrototypeCLITests.test_D5_emits_the_literal_identity_arithmetic_request) ... ok

----------------------------------------------------------------------
Ran 1 test in 0.622s

OK
```

## Tracer cycle 2: both positive-unit classes for `D=12`

The next literal test was RED because the identity-only path returned one
positive unit class instead of `[1,0]` and `[2,1]`:

```text
test_D12_emits_both_positive_unit_square_classes (__main__.PrototypeCLITests.test_D12_emits_both_positive_unit_square_classes) ... FAIL
test_D5_emits_the_literal_identity_arithmetic_request (__main__.PrototypeCLITests.test_D5_emits_the_literal_identity_arithmetic_request) ... ok

======================================================================
FAIL: test_D12_emits_both_positive_unit_square_classes (__main__.PrototypeCLITests.test_D12_emits_both_positive_unit_square_classes)
----------------------------------------------------------------------
AssertionError: Lists differ: [{'el[48 chars], 1]}] != [{'el[48 chars], 1]}, {'element': [2, 1], 'id': 'unit-1', 'si[14 chars] 1]}]

Second list contains 1 additional elements.
First extra element 1:
{'element': [2, 1], 'id': 'unit-1', 'signatures': [1, 1]}

----------------------------------------------------------------------
Ran 2 tests in 1.246s

FAILED (failures=1)
```

Enumerating unit exponent parities, filtering exact signatures, and
normalizing modulo the square of a free unit made that literal test GREEN.

## Tracer cycle 3: direct `ker(Sq)` eligibility for `D=40`

The next test required the canonical norm-two ideal and its positive square
generator.  The existing implementation returned only the identity class:

```text
test_D40_keeps_the_nontrivial_class_in_the_kernel_of_Sq (__main__.PrototypeCLITests.test_D40_keeps_the_nontrivial_class_in_the_kernel_of_Sq) ... FAIL

======================================================================
FAIL: test_D40_keeps_the_nontrivial_class_in_the_kernel_of_Sq (__main__.PrototypeCLITests.test_D40_keeps_the_nontrivial_class_in_the_kernel_of_Sq)
----------------------------------------------------------------------
AssertionError: Lists differ: [{'id[83 chars] 1]]}] != [{'id[83 chars] 1]]}, {'id': 'ideal-1', 'norm': 2, 'ordinary_[49 chars]1]]}]

Second list contains 1 additional elements.
First extra element 1:
{'id': 'ideal-1', 'norm': 2, 'ordinary_class_coordinates': [1], 'row_hnf': [[2, 0], [0, 1]]}

----------------------------------------------------------------------
Ran 3 tests in 1.896s

FAILED (failures=1)
```

Bounded-norm class enumeration and the direct proof-enabled test that
`a^2=(alpha)` for a totally positive `alpha` made the three-test suite GREEN:

```text
test_D12_emits_both_positive_unit_square_classes ... ok
test_D40_keeps_the_nontrivial_class_in_the_kernel_of_Sq ... ok
test_D5_emits_the_literal_identity_arithmetic_request ... ok

----------------------------------------------------------------------
Ran 3 tests in 1.889s

OK
```

## Tracer cycle 4: the Cartesian HM set for `D=60`

The four arithmetic records were already exact when the next test required
the public document to state the corrected ideal-class condition explicitly:

```text
test_D60_emits_the_literal_cartesian_HM_representatives (__main__.PrototypeCLITests.test_D60_emits_the_literal_cartesian_HM_representatives) ... ERROR

======================================================================
ERROR: test_D60_emits_the_literal_cartesian_HM_representatives (__main__.PrototypeCLITests.test_D60_emits_the_literal_cartesian_HM_representatives)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/edgarcosta/projects/genus2isogenies/docs/superpowers/rm-goal/test_prototype.py", line 152, in test_D60_emits_the_literal_cartesian_HM_representatives
    document["hm_structure"],
KeyError: 'hm_structure'

----------------------------------------------------------------------
Ran 4 tests in 2.496s

FAILED (errors=1)
```

Adding `hm_structure`, including the literal
`ker(Sq:Cl(F)->Cl+(F))` condition and the Cartesian count, made the tracer
GREEN.

## Acceptance hardening

The remaining public tests were added incrementally for the ordinary `C4`
case (`D=145`), the ineligible ordinary 2-class (`D=205`), the `C3`
image-versus-kernel trap (`D=229`), nonmaximal and invalid-discriminant
rejections, malformed CLI JSON, byte stability, mandatory run checks, and the
bounded genus-theory count oracle.  Their individual intermediate terminal
chunks were not retained, so they are not presented as exact RED quotations.

## Final GREEN

The orchestrator reran the complete suite from the repository root after the
implementation handoff:

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
Ran 13 tests in 22.853s

OK
```

Exit status was `0`.

## Post-review conformance cycle

The parallel Spec review found two partial contract mismatches: an HM record
linked to its ideal data only through `ideal_id`, and a nonmaximal rejection
reused success-like proof metadata.  Tests for the literal frozen behavior
were added first.  The full command was RED with exit `1`:

```text
test_D60_emits_the_literal_cartesian_HM_representatives ... ERROR
test_nonmaximal_order_is_a_validated_rejection_without_HM_data ... FAIL

KeyError: 'ideal'

AssertionError: 'EXACT_PROOF_ENABLED' != 'EXACT_ORDER_REJECTION_ONLY'

----------------------------------------------------------------------
Ran 13 tests in 23.054s

FAILED (failures=1, errors=1)
```

The implementation then embedded an exact copy of the referenced ideal record
inside each HM record and introduced a status-specific order-rejection theorem
boundary.  The same full command was GREEN:

```text
----------------------------------------------------------------------
Ran 13 tests in 23.110s

OK
```

The arithmetic enumeration itself was unchanged.  A separate clean agent and
the black-box/Magma validator reran after this remediation; their final
results are recorded in `RERUN_TRANSCRIPT.md` and `BRUTE_FORCE.md`.
