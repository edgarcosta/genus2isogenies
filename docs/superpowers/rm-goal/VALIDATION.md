# Validation ledger

## Verdict

The frozen maximal-order arithmetic prototype passes its acceptance suite, a
clean process rerun, a black-box exact reconstruction over 78 real quadratic
fields, a partial independent Magma oracle, and an adversarial review.  This
validates the emitted ideal/unit arithmetic.  It does **not** validate an
abelian-surface quotient, a polarized target, or the full RM isogeny graph.

Toolchain used:

```text
SageMath 10.8
Magma V2.29-8
latexmk 4.87
```

## Reproducible commands and observed results

All arithmetic commands are run from the repository root.

### Public acceptance suite

```sh
PYTHONHASHSEED=0 sage -python docs/superpowers/rm-goal/test_prototype.py -v
```

Final post-review independent clean-agent result:

```text
Ran 13 tests in 23.363s
OK
exit 0
```

The 13 tests use only the public process boundary.  They cover the seven
literal fields, both required rejection cases, malformed CLI input, all
mandatory run-level checks, byte stability, and every positive fundamental
discriminant through 80.  `RERUN_TRANSCRIPT.md` preserves the complete test
list, stream hashes, timing, and minimal environment.

### Direct successful request

```sh
PYTHONHASHSEED=0 sage -python docs/superpowers/rm-goal/prototype.sage \
  enumerate --field-discriminant 60 --order-conductor 1
```

Observed process result:

```text
exit=0
status=OK_ARITHMETIC_ONLY
ordinary_class_group=[2]
narrow_class_group=[2,2]
positive_unit_class_count=2
eligible_ordinary_class_count=2
record_count=4
alpha=[[1,0],[4,1],[2,0],[8,2]]
stdout_bytes=4921
stdout_sha256=80d0bb87f5c7e6851f4a262ae5de0fb6377a27e7c3f6472f8999fcdc514ab9ba
stderr_bytes=0
```

Each of the four records contains an exact local copy of its referenced ideal
ID, norm, row HNF, and ordinary class coordinates.  Two new processes produced
the same 4,921 bytes and the same SHA-256.

### Validated rejections

```sh
PYTHONHASHSEED=0 sage -python docs/superpowers/rm-goal/prototype.sage \
  enumerate --field-discriminant 5 --order-conductor 2
```

```text
exit=2
status=UNSUPPORTED_NONMAXIMAL_ORDER
field_discriminant=5
order_discriminant=20
index_in_maximal_order=2
arithmetic_result=EXACT_ORDER_REJECTION_ONLY
class_unit_enumeration=NOT_RUN
hm_records absent
```

```sh
PYTHONHASHSEED=0 sage -python docs/superpowers/rm-goal/prototype.sage \
  enumerate --field-discriminant 20 --order-conductor 1
```

```text
exit=2
status=INVALID_FIELD_DISCRIMINANT
hm_records absent
```

```sh
PYTHONHASHSEED=0 sage -python docs/superpowers/rm-goal/prototype.sage \
  enumerate --field-discriminant not-an-integer
```

```text
exit=2
status=INVALID_ARGUMENTS
stderr="invalid arguments: argument --field-discriminant: invalid int value: 'not-an-integer'"
```

### Bounded black-box oracle

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 sage -python \
  docs/superpowers/rm-goal/bruteforce_validation.sage \
  --max-discriminant 260 --require-magma
```

Observed exact summary:

```json
{"fundamental_discriminants_checked":78,"magma_cross_check":"PASS_7_DIAGNOSTIC_FIELDS","max_discriminant":260,"proof_mode":true,"prototype_access":"PUBLIC_CLI_ONLY","reconstructed_hm_records":148,"status":"PASS"}
```

Exit status was `0` and stderr was empty.  `BRUTE_FORCE.md` preserves the full
one-line JSON output, including every diagnostic field.

The validator never imports the prototype.  For every JSON record it creates
fresh exact Sage objects and checks the ideal HNF, norm, class coordinate,
integrality, invertibility, `a^2=(alpha)`, both signs, norm identity, unit
square class, and pair inequivalence.  It also computes `ker(Sq)` directly.

### Paper build

From `private/tex/_worktrees/rm-goal/`:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Observed result:

```text
exit=0
Output written on main.pdf (21 pages, 536027 bytes).
pdf_sha256=ede2e538aca4e8568f3bf428756e72900feb51be0bd96a8f146ad905c1967795
```

The only unresolved citations were the baseline keys `Gro72` and `BLR90`.
The added `GonzalezGuardiaRotger2005` entry resolved.  The complete retained
output is in `PAPER_LATEXMK.md`.  A final full orchestrator rerun also passed
with 21 pages and 536027 bytes.  Its PDF hash differed because of build
metadata; deterministic PDF bytes are not an acceptance criterion.

## Literal fixture matrix

`sign` is the size of the unit-signature image and `U+` is
`#(O_F^{x,+}/(O_F^x)^2)`.

| `D` | `Cl(F)` | `Cl+(F)` | sign | `#ker(Sq)` | `U+` | records | Exact purpose |
|---:|:---:|:---:|---:|---:|---:|---:|:---|
| 5 | 1 | 1 | 4 | 1 | 1 | 1 | negative-norm unit; identity only |
| 12 | 1 | C2 | 2 | 1 | 2 | 2 | unequal class groups; two polarization classes |
| 40 | C2 | C2 | 4 | 2 | 1 | 2 | admissible nontrivial ordinary 2-class |
| 60 | C2 | C2 x C2 | 2 | 2 | 2 | 4 | nontrivial ideal/unit Cartesian product |
| 145 | C4 | C4 | 4 | 2 | 1 | 2 | only the ordinary order-two subgroup survives |
| 205 | C2 | C4 | 2 | 1 | 2 | 2 | ordinary 2-class is principal-square but not positive-principal-square |
| 229 | C3 | C3 | 4 | 1 | 1 | 1 | squaring image is all C3, kernel is trivial |

For the seven fields, Magma independently reproduced both ordinary and narrow
class-group invariant factors.  The exact canonical ideals and `alpha`
values were reconstructed by the black-box validator but were not separately
recomputed in Magma.

## Independent-validation labels

| Claim | Evidence | Label |
|:---|:---|:---|
| HM record cardinality | Elementary trial-division value `2^(omega(D)-1)` for 78 fields | **INDEPENDENT FORMULATION** |
| Same cardinality | Fresh proof-enabled `#Cl+(F)[2]` and direct `#ker(Sq) * #U+/U^2` | **BLACK-BOX, SHARED SAGE/PARI ENGINE** |
| Ordinary/narrow invariant factors | Magma on `D=5,12,40,60,145,205,229` | **INDEPENDENT CAS, PARTIAL** |
| Exact HNF, `a^2=(alpha)`, positivity, norm, and pair inequivalence | Fresh objects in a separate validator process | **BLACK-BOX, SHARED SAGE/PARI ENGINE** |
| Literal canonical representatives | Seven fixed expected JSON fixtures plus exact reconstruction | **NOT INDEPENDENTLY VALIDATED BY A SECOND CAS** |
| Unit-square representatives | Exact unit coordinates and signatures in the black-box validator | **BLACK-BOX, SHARED SAGE/PARI ENGINE** |
| Deterministic serialization | Two fresh minimal-environment processes and a byte comparison | **INDEPENDENT PROCESS RERUN** |
| Exit-3 proof/internal failure variants | Static routing only; no safe fault injection was added | **NOT INDEPENDENTLY VALIDATED** |
| Fixed-surface polarization theorem | Primary-source audit of GGR Theorem 2.10 and the direct automorphism action | **PRIMARY-SOURCE VALIDATED** |
| LaTeX integrity | Required `latexmk` invocation | **BUILD VALIDATED, NOT MATHEMATICAL VALIDATION** |
| Quotient PPAS, transported RM, target equality, and full graph | No implementation or oracle in this slice | **NOT INDEPENDENTLY VALIDATED / OUT_OF_SCOPE** |

Agreement between local `note.tex` and arXiv:2405.19820v1 is explicitly not an
independent check; `ARXIV_DIFF.log` records only the comparison.

## Adversarial disposition

`ADVERSARIAL_REVIEW.md` gives every attack an explicit disposition.  No
blocker was found.  The decisive counterexamples are:

- `D=40`: an image-of-squaring rule would miss an eligible class;
- `D=205`: ordinary principality of the square does not imply a totally
  positive generator;
- `D=229`: a surjective squaring image would add two ineligible classes; and
- `D=5,f=2` versus `D=20,f=1`: an order discriminant is not a field
  discriminant.

Two prototype-only limitations remain visible rather than silently fixed:

1. `--help` is an argparse text metacommand, outside the JSON result grammar.
2. The 10,000-step normalization guard is operational, not a proved input
   bound; exceeding it produces an explicit exit-3 state rather than a false
   result.

Formal review also found two partial contract mismatches, both resolved before
the final reruns: HM records now include record-local ideal provenance, and a
nonmaximal result now uses status-specific theorem metadata.  The production
slice in `NEXT_SLICE.md` preserves those corrections and resolves the two
remaining limitations at the type and API boundary.
