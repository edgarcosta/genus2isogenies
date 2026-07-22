# Two-axis code review

Review checkpoint: `c29edc4e108b28aa5c2b66c9ed862c9409ff7f91`

Fixed point: `1c8d56c606443ba622bab8e096e78d3bf630b870`

Diff command:

```sh
git diff 1c8d56c606443ba622bab8e096e78d3bf630b870...c29edc4e108b28aa5c2b66c9ed862c9409ff7f91
```

Commit list at review launch:

```text
c29edc4 Add audited maximal RM arithmetic slice
```

The two axes were run by separate read-only subagents.  Standards sources were
the repository `AGENTS.md`, `README.md`, and the fixed Fowler-smell baseline.
The Spec sources were
`docs/superpowers/goal-rm-polarizations-and-isogenies.md` and
`FROZEN_SLICE.md`.

## Standards

One hard violation was reported:

- `test_prototype.py` declared `#!/usr/bin/env python3`, although repository
  `AGENTS.md` requires all Python scripts to run under `sage -python`.  The
  documented commands were correct and the file was non-executable, but the
  header still advertised a forbidden invocation.

Judgement-call smells:

- **Primitive Obsession:** closed statuses and theorem-boundary concepts are
  strings spread across JSON dictionaries and exception routing.  The typed
  result union in `NEXT_SLICE.md` is the intended production repair.
- **Divergent Change:** the throwaway `prototype.sage` combines arithmetic,
  canonicalization, serialization, parsing, and error transport.  This is
  accepted for the prototype; the production slice separates the API/result
  model from its CLI adapter.
- **Duplicated Code:** `_positive_unit_classes` and
  `_unit_parity_elements` repeat the real-quadratic unit-basis validation.
  This is a small refactor opportunity, not a correctness issue.
- **Duplicated Code, deliberately retained:** elementary factorization and
  fundamental-discriminant classification occur in both the public tests and
  black-box validator.  Extracting them would weaken the validator's stated
  independent-formulation role, so this duplication is justified.

No other documented-standard violation was found.  Sage commands, nested
paper-worktree location, unrelated-file preservation, and paper patch scope
conformed.

## Spec

No blocker or scope creep was reported.  Two nonblocking partial requirements
were found:

- Frozen contract: “Every record includes the ideal norm/HNF, ordinary class
  coordinates.”  Each `hm_records[]` item had only `ideal_id`; the promised
  values lived in `hm_ideal_classes[]`.  They were recoverable, but not
  literally record-local.
- Frozen contract: the result must expose an “arithmetic-only theorem
  boundary.”  `UNSUPPORTED_NONMAXIMAL_ORDER` inherited
  `arithmetic_result=EXACT_PROOF_ENABLED` even though class/unit enumeration
  was intentionally not run.  The rejection and absence of HM data were
  otherwise correct.

The reviewer confirmed all terminal artifacts, the isolated two-file paper
diff and successful build, the arithmetic-only scope, and the absence of any
forbidden graph/worktree changes.

## Remediation and rerun

All hard/partial findings were addressed test-first:

1. the plain-Python test shebang was removed;
2. every HM record now embeds exactly the ideal record named by `ideal_id`;
3. the nonmaximal theorem boundary now says
   `EXACT_ORDER_REJECTION_ONLY`, `class_unit_enumeration=NOT_RUN`, and
   `geometric_interpretation=NOT_EVALUATED`.

The pre-fix suite failed at the two new assertions; the post-fix suite passed
all 13 tests.  The separate clean agent then confirmed byte-stable final
output, and the separate black-box/Magma agent reconfirmed 78 fields, 148
records, and Magma `7/7`.  The adversarial agent checked all 14 literal fixture
records for exact local/global ideal equality.  The judgement-call prototype
smells remain documented for production rather than being refactored after
Wave 3.

Summary: Standards had 1 hard finding (resolved) and 4 judgement-call observations; Spec had 2 nonblocking partial findings (resolved), with no blocker on either axis.
