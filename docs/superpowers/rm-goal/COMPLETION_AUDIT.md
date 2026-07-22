# Terminal requirement audit

Date: 2026-07-22 UTC

Verdict: **PASS**.  Every terminal artifact and scoped implementation
requirement is present.  The result is the maximal-order exact arithmetic
slice plus one fixed-surface paper proposition, not a full RM graph.

## Scope and ordering

- `BASELINE.md` was written before every other RM-goal artifact and records
  root, `private/`, and current `private/tex` status and commits.
- Positive claims use a geometrically simple surface over `k`, a chosen
  `iota:O_F -> End_k(A)` with image exactly `End_k(A)`, maximal real-quadratic
  `O_F`, and a reference principal polarization over `k`.
- `End_k`, `End_kbar`, the chosen embedding, and the Rosati-fixed subring are
  distinguished throughout.
- Nonmaximal orders, RM not defined over `k`, larger rational endomorphism
  rings, products, QM, matrix algebras, and geometric extra-polarization
  classification remain typed future obligations or out of scope.
- `FROZEN_SLICE.md` was completed after the three Wave 1 audits and before
  prototype implementation.  The implementation did not expand to a quotient
  engine, graph traversal, `Q`, `L_1`, or `L_2`.

## Sources and mathematical audit

- `private/tex/note.tex` was read completely; the specified RM,
  fundamental-isogeny, complex-building-block, and algorithm sections of the
  paper plus both planning files were audited.
- The typical implementation was used only for interface/oracle comparison.
- `ARXIV_DIFF.log` records the one comparison with arXiv:2405.19820v1 and
  correctly treats it as non-independent.  It found editorial/bibliographic,
  not mathematical, differences.
- `WAVE1_MATH_AUDIT.md`, `WAVE1_SOURCE_AUDIT.md`, and
  `WAVE1_IMPLEMENTATION_INVENTORY.md` have nonoverlapping Wave 1 ownership.
- Accepted theorems are tied to checked primary sources: Mumford for
  polarization/descent, González--Guàrdia--Rotger Theorem 2.10 for the
  fixed-surface classification, van der Geer for `Sq`/Hurwitz--Maass data,
  Hauffe-Waschbüsch--Krieg for the classical index, and Faltings only for the
  finiteness statement it actually proves.
- The frozen correction is
  `ker(Sq:Cl(F)->Cl+(F))`, extended by `O_F^{x,+}/(O_F^x)^2`; it does not use
  the note's incorrect image-of-squaring phrase.

## Wave 2 implementation and paper

- `prototype.sage` accepts exactly a fundamental field discriminant and order
  conductor, rejects a nonmaximal order, enables proofs globally, emits
  deterministic canonical JSON, checks every ideal square/generator/sign/norm,
  and keeps geometric degree labels theorem-conditional.
- Every final HM record contains local ideal ID, norm, row HNF, and ordinary
  coordinates equal to its referenced canonical ideal.
- `test_prototype.py` contains 13 public-process tests.  The TDD record includes
  actual RED/GREEN output and the post-review conformance cycle.
- The paper edit lives only in the fresh nested worktree
  `private/tex/_worktrees/rm-goal/`, whose parent repository locally excludes
  `/_worktrees/`.
- The paper branch changes only `main.tex` and `genus2classes.bib`; it adds one
  scoped proposition and one primary bibliography entry.  No HM/full-graph
  theorem was inserted.
- `PAPER_DIFF.patch` is byte-equal to the baseline-to-paper-branch diff and
  passes `git apply --check` in the clean current paper checkout.
- The required `latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex`
  command passed twice, including a final full rebuild.  Only the baseline
  `Gro72`/`BLR90` warnings remain.

## Wave 3 and post-review validation

- Separate clean-rerun, brute-force, and adversarial subagents produced
  `RERUN_TRANSCRIPT.md`, `BRUTE_FORCE.md` plus its validator, and
  `ADVERSARIAL_REVIEW.md`.
- After formal review remediation, those same independent owners reran on the
  amended files rather than relying on stale evidence.
- Final clean suite: 13/13, exit `0`; two `D=60` outputs were byte-identical at
  4,921 bytes with SHA-256
  `80d0bb87f5c7e6851f4a262ae5de0fb6377a27e7c3f6472f8999fcdc514ab9ba`.
- Final black-box oracle: 78 fundamental fields through `D=260`, 148 exact
  reconstructed records, elementary genus count PASS, direct narrow
  2-torsion/count PASS, and Magma ordinary/narrow invariants PASS on all seven
  diagnostic fields.
- The adversarial review attacks unit duplicates, omitted/ineligible ideal
  classes, nonprincipal ideals, positive-generator signs, order input,
  embedding dependence, and `A` versus PPAS/target identity.  Every finding
  has a disposition and no blocker remains.
- Evidence sharing Sage/Pari, exit-3 routing, and every geometric/analytic
  claim are labelled with their actual independence limits in `VALIDATION.md`.

## Formal review and next slice

- `CODE_REVIEW.md` preserves separate Standards and Spec reports.  The one
  hard Standards issue and both partial Spec issues were fixed test-first.
- `REPORT.md` answers every required input/data-layer, finite-orbit,
  polarization, quotient-return, phase-boundary, algorithm-step, and
  future-interface question.
- `NEXT_SLICE.md` specifies one API,
  `enumerate_maximal_rm_arithmetic(field_discriminant, order_conductor)`, with
  destination files, immutable result variants, nine fixtures, failure states,
  five tiny commits, and executable acceptance criteria.  It explicitly makes
  no change to `genus2isogenies.py`.

## Changed paths and repository state

Root commits:

```text
c29edc4 Add audited maximal RM arithmetic slice
4ee03c5 Address maximal RM review findings
```

Every root path changed since baseline
`1c8d56c606443ba622bab8e096e78d3bf630b870` is under
`docs/superpowers/rm-goal/`.  All unrelated untracked paths recorded in the
baseline remain untouched.

Paper commit on local branch `rm-goal-paper`:

```text
c6f168d Document fixed-surface RM polarizations
```

The current `private/tex` checkout remains clean at its original commit.  The
pre-existing modification visible in `_worktrees/gluing` during the terminal
audit was not opened or changed by this work; neither it nor
`_worktrees/isogeny-primes-engine` appears in either goal diff.

No push, PR, or full RM graph implementation was performed.
