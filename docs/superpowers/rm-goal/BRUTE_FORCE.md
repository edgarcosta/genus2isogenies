# Independent bounded arithmetic validation

## Result

PASS.  Every positive fundamental real-quadratic discriminant
`2 <= D <= 260` was checked: 78 fields and 148 reported HM records.  For each
field, the number of `hm_records` agreed both with the elementary genus-theory
count

```text
2^(omega(D)-1)
```

and with the proof-enabled Sage computation of `#Cl+(F)[2]`.  All reported
ideals, elements, and equivalence-class identifiers survived fresh exact
reconstruction.  Magma independently matched the ordinary and narrow class
group invariant factors in all seven diagnostic fields.

This is an arithmetic validation only.  It says nothing about construction or
uniqueness of quotient PPAS targets.

Post-review rerun (2026-07-22): PASS again with the exact command below after
the CLI schema was strengthened to embed the complete ideal payload in every
HM record.  The amended validator checked all 148 record-local payloads for
exact equality with the canonical ideal record referenced by `ideal_id`.

## Reproduction

Run from the repository root:

```sh
PYTHONHASHSEED=0 sage -python docs/superpowers/rm-goal/bruteforce_validation.sage \
  --max-discriminant 260 --require-magma
```

Environment used:

```text
SageMath version 10.8, Release Date: 2025-12-18
Magma V2.29-8
```

Exact stdout from the successful run (one JSON line):

```json
{"diagnostic_fields":{"12":{"class_invariants":[],"direct_sq_image_size":1,"direct_sq_kernel_size":1,"hm_record_count":2,"narrow_class_invariants":[2],"narrow_two_torsion_count":2,"positive_unit_class_count":2,"signature_image_size":2},"145":{"class_invariants":[4],"direct_sq_image_size":2,"direct_sq_kernel_size":2,"hm_record_count":2,"narrow_class_invariants":[4],"narrow_two_torsion_count":2,"positive_unit_class_count":1,"signature_image_size":4},"205":{"class_invariants":[2],"direct_sq_image_size":2,"direct_sq_kernel_size":1,"hm_record_count":2,"narrow_class_invariants":[4],"narrow_two_torsion_count":2,"positive_unit_class_count":2,"signature_image_size":2},"229":{"class_invariants":[3],"direct_sq_image_size":3,"direct_sq_kernel_size":1,"hm_record_count":1,"narrow_class_invariants":[3],"narrow_two_torsion_count":1,"positive_unit_class_count":1,"signature_image_size":4},"40":{"class_invariants":[2],"direct_sq_image_size":1,"direct_sq_kernel_size":2,"hm_record_count":2,"narrow_class_invariants":[2],"narrow_two_torsion_count":2,"positive_unit_class_count":1,"signature_image_size":4},"5":{"class_invariants":[],"direct_sq_image_size":1,"direct_sq_kernel_size":1,"hm_record_count":1,"narrow_class_invariants":[],"narrow_two_torsion_count":1,"positive_unit_class_count":1,"signature_image_size":4},"60":{"class_invariants":[2],"direct_sq_image_size":1,"direct_sq_kernel_size":2,"hm_record_count":4,"narrow_class_invariants":[2,2],"narrow_two_torsion_count":4,"positive_unit_class_count":2,"signature_image_size":2}},"fundamental_discriminants_checked":78,"magma_cross_check":"PASS_7_DIAGNOSTIC_FIELDS","max_discriminant":260,"proof_mode":true,"prototype_access":"PUBLIC_CLI_ONLY","reconstructed_hm_records":148,"status":"PASS"}
```

Exit status was `0`; stderr was empty.

## What the validator checks

The validator treats `prototype.sage` as a black box.  It starts a new process
through the frozen command-line interface for every field and never imports
prototype code.  It then constructs a fresh number field, maximal order, class
group, narrow class group, unit group, and real embeddings with Sage proof mode
enabled.

For every returned document it checks:

- the canonical polynomial, field discriminant, maximal-order claim, and both
  class-group invariant lists;
- the unit signature image and the independently derived cardinality
  `#(U+/U^2) = 4 / #image(sign)`;
- that every reported positive-unit representative is a unit, is totally
  positive, and is distinct modulo unit squares;
- that every ideal row-HNF reconstructs an integral invertible ideal, its
  determinant equals its declared norm, and its ordinary class coordinates are
  correct and nonrepeated;
- the kernel and image sizes of `Sq: Cl(F) -> Cl+(F)`, computed directly from
  ideal quotients and the signs of proof-certified principal generators rather
  than by trusting the JSON diagnostics;
- for every record, the exact equality `a^2 = (alpha)`, total positivity of
  `alpha`, `Norm(alpha) = Norm(a)^2`, and the declared graph-weight/degree norm
  identities;
- for every record, exact JSON equality between its record-local `ideal`
  payload and the canonical top-level ideal record referenced by `ideal_id`;
- uniqueness of unit, ideal, and record IDs, uniqueness of every
  `(ideal_id, unit_class_id)`, and pairwise inequivalence under
  `(a,alpha) ~ (x*a,x^2*alpha)` using exact principal-ideal and unit-square
  tests; and
- the three independent cardinality expressions
  `#records = 2^(omega(D)-1) = #Cl+(F)[2]
  = #ker(Sq) * #(U+/U^2)`.

The field-range classifier and the first count use elementary squarefreeness
and trial-division routines implemented in the validator, not Sage's
fundamental-discriminant or factorization routines.  The selected range is
also checked against Sage's classifier as a consistency assertion.

## Required diagnostic coverage

| D | `Cl(F)` | `Cl+(F)` | sign-image size | `#ker(Sq)` | `#im(Sq)` | `#U+/U^2` | records | phenomenon |
|---:|:---:|:---:|---:|---:|---:|---:|---:|:---|
| 5 | 1 | 1 | 4 | 1 | 1 | 1 | 1 | a norm-`-1` unit gives the full sign image |
| 12 | 1 | C2 | 2 | 1 | 1 | 2 | 2 | no norm-`-1` unit; ordinary and narrow class numbers differ |
| 40 | C2 | C2 | 4 | 2 | 1 | 1 | 2 | the nontrivial ordinary 2-class is admissible |
| 60 | C2 | C2 x C2 | 2 | 2 | 1 | 2 | 4 | unequal class numbers and a nontrivial ideal/unit product |
| 145 | C4 | C4 | 4 | 2 | 2 | 1 | 2 | ordinary C4; only its order-two subgroup lies in the kernel |
| 205 | C2 | C4 | 2 | 1 | 2 | 2 | 2 | the nontrivial ordinary 2-class is ineligible |
| 229 | C3 | C3 | 4 | 1 | 3 | 1 | 1 | squaring is surjective but its kernel is trivial |

Thus both unit-signature regimes, unequal ordinary/narrow class groups,
admissible and inadmissible ordinary 2-classes, ordinary C4, and the C3
image-versus-kernel trap are exercised explicitly.

## Independence boundary

There are three levels of independence:

1. The genus-theory value uses an elementary factor counter that shares no
   implementation with the prototype.
2. The certificate checks use a separate process and fresh arithmetic objects,
   and they do not import prototype internals.  They nevertheless use the same
   installed Sage 10.8/Pari arithmetic engine as the prototype.  Consequently,
   this layer is an independent implementation and black-box check, but not an
   independent CAS oracle.
3. Installed Magma V2.29-8 independently recomputed `Cl(F)` and `Cl+(F)`
   invariant factors for exactly `D = 5, 12, 40, 60, 145, 205, 229`.  Magma did
   **not** inspect the prototype JSON, reconstruct its chosen ideals or
   `alpha` values, test pair equivalence, compute unit-signature representatives,
   or certify the analytic/geometric interpretation.  Those claims remain in
   levels 1--2 only.

No result here validates target-surface construction, polarization transport,
descent/twist data, or duplicate-free PPAS targets.
