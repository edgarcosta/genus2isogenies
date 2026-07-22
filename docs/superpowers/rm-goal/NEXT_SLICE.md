# Next production slice: exact maximal-RM arithmetic module

## Unequivocal recommendation

Extract only the proof-enabled maximal-order arithmetic enumerator into a new
production module.  Do not connect it to graph traversal or to the existing
modular-equation wrappers in this slice.

The single public Python API is:

```python
def enumerate_maximal_rm_arithmetic(
    field_discriminant: int,
    order_conductor: int = 1,
) -> RMArithmeticOutcome:
    ...
```

It accepts arithmetic field/order data only.  It does not accept a curve,
surface, polarization, embedding, period matrix, or prime list.  The return is
a closed tagged union, not a dictionary whose fields change silently by
status.

## Destination files

- `rm_arithmetic.py`: immutable input/result types, exact enumeration,
  invariants, and canonical JSON serialization.
- `tests/test_rm_arithmetic.py`: public-API and serialization tests run under
  `sage -python`.
- `tests/fixtures/rm_arithmetic/`: committed golden JSON for the seven
  acceptance fields and the two rejection cases.

Do not edit `genus2isogenies.py` yet.  Its current curve-only vertex loses
both the polarization and the RM embedding, so an adapter there would create
an unsafe seam.  The pre-existing untracked `tests/` contents recorded in
`BASELINE.md` must be inspected for ownership before adding only the named
paths; no other existing file in that directory is part of this slice.

## Types and invariants

Use immutable dataclasses (or equivalently explicit immutable records):

```text
RMArithmeticSuccess
  status = OK_ARITHMETIC_ONLY
  field: CanonicalRealQuadraticField
  order: MaximalOrderCertificate
  ordinary_class_group: ClassGroupDiagnostics
  narrow_class_group: ClassGroupDiagnostics
  unit_classes: tuple[PositiveUnitClass, ...]
  ideal_classes: tuple[HMIdealClass, ...]
  pairs: tuple[HMPair, ...]
  checks: ArithmeticChecks

RMInputRejected
  status: INVALID_FIELD_DISCRIMINANT
        | INVALID_ORDER_CONDUCTOR
        | UNSUPPORTED_NONMAXIMAL_ORDER
  input: RMOrderInput
  order_diagnostics: optional[OrderDiagnostics]

RMArithmeticUnknown
  status: PROOF_FAILURE_CLASS_GROUP
        | UNKNOWN_POSITIVE_GENERATOR
        | INTERNAL_ARITHMETIC_INCONSISTENCY
  diagnostic: str
```

Only `RMArithmeticSuccess` may carry unit, ideal, pair, or theorem-conditional
degree data.  This removes the prototype's overbroad theorem payload on a
nonmaximal rejection.  Help belongs to a separate CLI adapter and is
explicitly exempt from the result-JSON contract.

Every success must retain these invariants:

- canonical field basis `[1,w]` and maximal-order equality;
- proof-enabled class, narrow-class, unit, and principal-ideal computations;
- least `(norm,row-HNF)` integral representative within the real-quadratic
  Minkowski bound for every ordinary class;
- direct eligibility test `a^2=(alpha)` with `alpha>>0`;
- exact `Norm(alpha)=Norm(a)^2` and integral/invertible checks;
- representatives of `O_F^{x,+}/(O_F^x)^2` distinct in unit exponent parity;
- pair representatives distinct under
  `(a,alpha)~(beta*a,beta^2*alpha)`; and
- stable ordering and canonical JSON independent of Python hash order.

Replace the prototype's iterative 10,000-step normalization with a reduction
whose exponent is selected directly from certified real-algebraic bounds, or
return a named resource/precision-unknown result.  Do not encode a heuristic
iteration ceiling as an arithmetic inconsistency.

## Failure and excluded states

The production arithmetic union must distinguish exactly:

- `INVALID_FIELD_DISCRIMINANT`;
- `INVALID_ORDER_CONDUCTOR`;
- `UNSUPPORTED_NONMAXIMAL_ORDER`;
- `PROOF_FAILURE_CLASS_GROUP` (including narrow-class/unit certification);
- `UNKNOWN_POSITIVE_GENERATOR`; and
- `INTERNAL_ARITHMETIC_INCONSISTENCY`.

A later geometric adapter, in a different module and commit series, must add
distinct states for:

- `UNSUPPORTED_RM_NOT_DEFINED_OVER_BASE`;
- `UNSUPPORTED_LARGER_RATIONAL_ENDOMORPHISM_RING`;
- `UNKNOWN_GEOMETRIC_ENDOMORPHISMS`;
- `UNKNOWN_ANALYTIC_QUOTIENT`; and
- `UNKNOWN_PPAS_COLLISION`.

Those future states require a marked input `(A,lambda,iota)`.  They must never
be inferred from the discriminant-only arithmetic API.

## Fixed fixtures

The production tests must preserve literal expectations for:

| Input | Required phenomenon |
|:---|:---|
| `D=5,f=1` | one unit, ideal, and pair |
| `D=12,f=1` | two positive-unit square classes |
| `D=40,f=1` | eligible nontrivial ordinary 2-class |
| `D=60,f=1` | two ideals times two units gives four pairs |
| `D=145,f=1` | only the order-two subgroup of ordinary C4 |
| `D=205,f=1` | nontrivial ordinary 2-class rejected by signs |
| `D=229,f=1` | trivial kernel despite surjective squaring on C3 |
| `D=5,f=2` | typed nonmaximal rejection, no HM payload |
| `D=20,f=1` | typed invalid-field rejection, no HM payload |

Retain the black-box range check through `D=260` as a slower acceptance test,
including the elementary `2^(omega(D)-1)` oracle.  Magma is a validation job,
not a runtime dependency of `rm_arithmetic.py`.

## Tiny commit sequence

1. **Add result types and validation gates.**  Implement only canonical field
   construction, the maximal-order gate, and the three input rejection
   variants; add `D=5,f=2` and `D=20,f=1` tests.
2. **Add the positive-unit quotient.**  Implement exact sign/parity classes
   and canonical serialization; add the `D=5` and `D=12` goldens.
3. **Add direct `ker(Sq)` representatives.**  Add bounded ordinary-class
   enumeration, positive square generators, pair records, and literal
   `D=40,60,145,205,229` goldens.
4. **Add invariant and black-box acceptance tests.**  Reconstruct every pair,
   check byte stability and the `D<=260` genus count, and run the seven-field
   Magma oracle as an optional validation command.
5. **Document the geometric adapter boundary.**  Add no graph code; record the
   future marked-PPAS request/response protocol and excluded states.

Each commit is independently green under the tests introduced with it.  No
commit changes `genus2isogenies.py` or claims a computed quotient.

## Executable acceptance criteria

From the repository root, all of the following must pass:

```sh
PYTHONHASHSEED=0 sage -python -m unittest tests/test_rm_arithmetic.py -v
```

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 sage -python \
  docs/superpowers/rm-goal/bruteforce_validation.sage \
  --max-discriminant 260 --require-magma
```

The unit suite must assert:

- exact equality with all nine golden outcomes;
- one-line, byte-stable canonical JSON for every outcome variant;
- no success-only fields on either rejection or unknown results;
- all mandatory exact checks on every success;
- direct pairwise inequivalence of emitted pairs; and
- no call to an analytic quotient, modular equation, curve reconstruction, or
  graph traversal function.

The bounded validator must report exactly:

```text
status=PASS
fundamental_discriminants_checked=78
reconstructed_hm_records=148
magma_cross_check=PASS_7_DIAGNOSTIC_FIELDS
```

Only after this slice lands should a separate design task specify the marked
`(A,lambda,iota)` quotient adapter and its exact polarized-isomorphism key.
