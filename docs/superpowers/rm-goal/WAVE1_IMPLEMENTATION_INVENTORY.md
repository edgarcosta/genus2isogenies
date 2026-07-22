# Wave 1 implementation inventory: maximal RM arithmetic

## Decision

The smallest safe prototype is an **arithmetic candidate enumerator**, not an
isogeny-graph implementation.  It accepts a real quadratic order, rejects it
unless it is maximal, and returns canonical exact representatives of

1. totally positive unit classes
   \(\mathcal O_F^{\times,+}/(\mathcal O_F^\times)^2\); and
2. Hurwitz--Maass pairs \((\mathfrak p,\alpha)\), modulo
   \((\mathfrak p,\alpha)\sim(\beta\mathfrak p,\beta^2\alpha)\), satisfying
   \(\mathfrak p^2=(\alpha)\) with \(\alpha\gg0\).

This is exact number-field arithmetic.  Interpreting a returned pair as a
rational polarized quotient is **conditional on the theorem audit**.  The
prototype must not claim to enumerate quotient surfaces, the full RM isogeny
class, or even distinct unmarked PPAS targets.

The scope assumed by any later geometric adapter is:

- `A/k` is geometrically simple;
- a particular embedding `iota: O_F -> End_k(A)` is part of the data;
- `O_F` is the maximal order of a real quadratic field;
- `End_k(A) = iota(O_F)`, not merely `O_F subset End_k(A)`;
- the chosen principal polarization `lambda` and its Rosati involution are
  explicit, and `iota(O_F)` is Rosati-fixed.

`End_k(A)`, `End_kbar(A)`, the chosen embedding `iota`, and
`End_k(A)^dagger_lambda` are separate fields in the future interface.  Known
larger geometric endomorphism rings, unknown RM descent, nonmaximal orders,
products, QM, and `M_2(Q)` are not silently accepted.

## Why the direct ideal test is the executable definition

The note first requires

\[
  \mathfrak p^2=(\alpha),\qquad \alpha\gg0
  \tag{HM}
\]

(`note.tex`, lines 228--243), but later asks for “classes in
`Cl^+(F)` that are squares of elements of `Cl(F)`” (lines 245--255).  The
latter phrase is not well-typed without a chosen map or section and is not an
implementation rule.

Define instead

\[
 H_F=\{c\in\mathrm{Cl}(F):\text{one/every integral ideal }\mathfrak p\in c
       \text{ satisfies (HM)}\}.
\]

The property is independent of the representative: replacing `p` by
`beta*p` replaces `alpha` by `beta^2*alpha`, and `beta^2` is totally positive.
The prototype computes `H_F` directly.  It does **not** select either
`Cl(F)[2]` or `Cl^+(F)[2]`.  Conditionally on the standard narrow-class exact
sequence, `H_F` is the image of `Cl^+(F)[2] -> Cl(F)`, but the code does not
need that reformulation.

This matters in `F` of discriminant 205:

- `Cl(F) = C2` and `Cl^+(F) = C4`;
- the nontrivial ordinary class squares to a principal ideal whose generator
  has mixed signs, and no unit changes it to a totally positive generator;
- hence `#H_F = 1`, although both ordinary and narrow 2-torsion have size 2.

For each canonical `p` in `H_F`, the positive generators of `p^2` modulo
unit squares form a torsor for
`O_F^(x,+)/(O_F^x)^2`.  Thus the arithmetic records are the finite product

\[
  H_F\times \mathcal O_F^{\times,+}/(\mathcal O_F^\times)^2.
\]

This has no duplicate **arithmetic pair classes**.  Two resulting quotients
could nevertheless be isomorphic after forgetting the source isogeny or the
RM marking; only the quotient engine can detect that, so it must return an
exact isomorphism/deduplication certificate.

## Exact Sage and CLI contract

The proposed file is `docs/superpowers/rm-goal/prototype.sage`, written as
ordinary Python importing `sage.all`, and run by:

```sh
PYTHONHASHSEED=0 sage -python docs/superpowers/rm-goal/prototype.sage \
  enumerate --field-discriminant 60 --order-conductor 1
```

No `proof=False` path is allowed.  Stdout is one deterministic JSON document;
diagnostics go to stderr; success is exit 0, rejected scope is exit 2, and an
unproved/internal failure is exit 3.  The command writes no files.

Input:

```text
field_discriminant : positive fundamental discriminant D
order_conductor    : positive integer f (default 1)
```

Use the canonical integral generator `w`:

```text
D odd:  w^2 - w + (1-D)/4 = 0
D even: w^2 - D/4 = 0
```

Hence `O_F = Z[1,w]` and the supplied order is `Z[1,f*w]`.  The exact in-memory
objects are a Sage `NumberField`, `Order`, `NumberFieldFractionalIdeal`, and
`NumberFieldElement`.  The gate `O == K.maximal_order()` occurs before class
or unit enumeration.

JSON uses only integers and arrays:

```text
field.defining_polynomial       ascending integer coefficients
field.discriminant              D
order.basis                     [[1,0],[0,1]] in [1,w]
ordinary_class_group.invariants invariant factors
narrow_class_group.invariants   diagnostic invariant factors only
unit_signature_image            sorted sign pairs
polarization_unit_classes[]     {id, element:[a,b], signatures}
hm_ideal_classes[]              {id, norm, row_hnf, ordinary_class_coordinates}
hm_records[]                    {id, ideal_id, unit_class_id, alpha:[a,b], checks}
```

`[a,b]` means `a+b*w`.  `row_hnf` is the row Hermite normal form of the ideal
lattice in the fixed basis `[1,w]`.  Ordinary class coordinates are diagnostic
only; IDs and sorting use `(norm,row_hnf)`.

Every output also contains:

```json
{
  "schema": "rm-hm-arithmetic-v1",
  "status": "OK_ARITHMETIC_ONLY",
  "surface_scope": "NOT_EVALUATED",
  "theorem_status": "GEOMETRIC_INTERPRETATION_CONDITIONAL"
}
```

## Deterministic enumeration

1. Validate `D`, construct the canonical polynomial and `O`, and reject
   `f != 1` (also verify equality with the maximal order rather than trusting
   the integer).
2. Compute ordinary and narrow class groups and the unit group with proofs.
   The narrow group is reported as a diagnostic, not used to select ideals.
3. Let `B=floor(sqrt(D)/2)`, the real-quadratic Minkowski bound.  Enumerate
   all integral ideals of norm at most `B`.  For each ordinary class choose the
   least `(norm,row_hnf)`.  Assert that every ordinary class was represented.
4. For each chosen ideal `p`, compute `p^2`.  Discard it unless it is principal.
   If a generator is `g`, enumerate the four exponent-parity classes in
   `O_F^x/(O_F^x)^2`; retain exactly those `u*g` that are totally positive.
   The ideal class belongs to `H_F` iff this list is nonempty.
5. Independently retain totally positive representatives among the same four
   unit parity classes to obtain `O_F^(x,+)/(O_F^x)^2`.
6. Canonicalize a positive generator modulo square units.  If `e` is a free
   unit generator, orient `eta=e^2` so its value at the larger real embedding
   is larger.  Multiply by a unique power of `eta` so that
   `sigma_+(a) >= sigma_-(a)` but
   `sigma_+(a/eta) < sigma_-(a/eta)`.  Use exact algebraic-real comparisons,
   not floating signs.  Sort coefficient pairs and remove exact duplicates.
7. Emit the Cartesian records in `(ideal norm, HNF, unit coefficients)` order.
   Include `(O_F,1)` and label it `identity_data`; do not confuse omission of
   the identity with a group computation.

The canonicalization is relative only to the fixed field presentation and
integral basis, not Sage's choice of class-group generators.  Conjugate ideals
are **not** identified: the chosen RM embedding makes them different marked
requests.

## Mandatory checks on every run

- the polynomial, maximal order, and order discriminants equal `D`;
- the Minkowski enumeration covers exactly `#Cl(F)` ordinary classes;
- all selected ideals are nonzero integral invertible ideals;
- for every record, exact ideal equality `p^2 == (alpha)`;
- `alpha` is integral and totally positive at both exact real embeddings;
- `Norm(alpha) == Norm(p)^2`;
- the expected quotient degree is `Norm(p)^2` and graph weight is `Norm(p)`
  (these two geometric interpretations remain theorem-conditional);
- unit representatives are totally positive and pairwise distinct modulo
  squares, checked in unit-group exponent coordinates;
- HM pair records are pairwise distinct under the stated ideal/unit
  equivalence;
- repeated invocations have byte-identical stdout.

## Explicit states

| Status | Meaning |
|---|---|
| `OK_ARITHMETIC_ONLY` | Exact ideal/unit data; no surface was inspected. |
| `INVALID_FIELD_DISCRIMINANT` | `D` is not a positive fundamental real-quadratic discriminant. |
| `UNSUPPORTED_NONMAXIMAL_ORDER` | The supplied order differs from `O_F`; include its conductor/index/discriminant. |
| `PROOF_FAILURE_CLASS_GROUP` | Sage could not certify a class, narrow-class, or unit computation. |
| `UNKNOWN_POSITIVE_GENERATOR` | A principal-square/sign computation did not certify either inclusion or exclusion. |
| `INTERNAL_ARITHMETIC_INCONSISTENCY` | Any mandatory equality/count check failed. |
| `UNSUPPORTED_SURFACE_SCOPE` | Future adapter: not geometrically simple, RM not defined over `k`, `End_k(A) != O_F`, or Rosati compatibility fails. |
| `UNKNOWN_GEOMETRIC_ENDOMORPHISMS` | Future adapter cannot distinguish the scoped RM case from a larger geometric endomorphism case. |
| `UNKNOWN_ANALYTIC_QUOTIENT` | Interval/recognition computation did not certify a quotient. |
| `UNKNOWN_PPAS_COLLISION` | A moduli collision lacks an exact polarized-isomorphism certificate. |

An exception traceback is never a mathematical negative answer.  In
particular, nonmaximal input is not coerced to the maximal order.

## Fixture matrix

These expectations were checked locally with Sage 10.8 using proof-enabled
class, narrow-class, and unit groups.

| `D` | `Cl` | `Cl+` | unit signature image size | `#H_F` | `#U+/U^2` | HM records | Purpose |
|---:|---|---|---:|---:|---:|---:|---|
| 5 | `1` | `1` | 4 | 1 | 1 | 1 | negative-norm unit; trivial action |
| 12 | `1` | `C2` | 2 | 1 | 2 | 2 | unequal class/narrow class; two polarizations |
| 40 | `C2` | `C2` | 4 | 2 | 1 | 2 | admissible nontrivial ideal 2-torsion |
| 60 | `C2` | `C2 x C2` | 2 | 2 | 2 | 4 | unequal groups plus admissible 2-torsion |
| 145 | `C4` | `C4` | 4 | 2 | 1 | 2 | only the order-two ordinary class survives |
| 205 | `C2` | `C4` | 2 | 1 | 2 | 2 | nontrivial ordinary 2-torsion fails (HM) |

Required rejection fixture:

```sh
PYTHONHASHSEED=0 sage -python docs/superpowers/rm-goal/prototype.sage \
  enumerate --field-discriminant 5 --order-conductor 2
```

It must exit 2 with `UNSUPPORTED_NONMAXIMAL_ORDER` (order discriminant 20),
and emit no HM data.  An invalid-discriminant fixture (`D=20,f=1`) separately
checks that an order discriminant is not accepted as a field discriminant.

## Principal polarizations on a fixed surface (conditional contract)

For fixed `(A,lambda,iota)` satisfying the scope, the intended correspondence
is

\[
  \epsilon\in O_F^{\times,+}\longmapsto
  \lambda_\epsilon=\lambda\circ\iota(\epsilon).
\]

Pullback by `u in Aut_k(A)=O_F^x` sends `epsilon` to
`u^dagger epsilon u = u^2 epsilon`, so the expected polarized-isomorphism
classes are `O_F^(x,+)/(O_F^x)^2`.  This statement must remain labelled
`CONDITIONAL_ON_POLARIZATION_CORRESPONDENCE` until the primary-source audit
certifies the precise field-of-definition and equivalence hypotheses.

Changing the reference polarization to `lambda_beta` changes Rosati by

\[
 x^{\dagger_{\lambda_\beta}}=\beta^{-1}x^{\dagger_\lambda}\beta.
\]

It remains the identity on the commutative ring `iota(O_F)` in this scope, but
it need not remain unchanged on a larger geometric endomorphism algebra.  The
prototype therefore reports unit classes relative to the supplied `lambda`;
it never treats an unpolarized abelian variety as a PPAS vertex.

## Contract required from an analytic quotient engine

Input is a marked source
`(A,lambda,iota,period_data)` plus an exact `(p,alpha)` record.  The engine must
return either an explicit status above or all of:

1. a certified kernel equal to `A[p]`, its Galois stability, and a quotient
   isogeny `phi` of the expected degree;
2. an exact PPAS `(B,lambda_B)` and a certificate
   `phi^vee lambda_B phi = lambda iota(alpha)`;
3. the transported embedding `iota_B` satisfying
   `phi iota(x) = iota_B(x) phi`, together with the target endomorphism-order
   status (continue only if it remains in scope);
4. certified rational modular invariants, an exact model over `k`, and twist
   resolution; a unit record with `p=O_F` must be handled as a
   polarization-only transition on the same underlying abelian variety;
5. a normalized period matrix, polarization form, RM analytic action, and
   precision certificate sufficient to launch the next request;
6. a canonical PPAS fingerprint plus an exact polarized-isomorphism
   certificate for every deduplication, explicitly recording whether the RM
   marking is preserved, conjugated, or forgotten;
7. trace/hash and degree data usable by the existing validation oracles.

Approximate Igusa invariants alone are insufficient.  Failure to recognize a
rational model is `UNKNOWN_ANALYTIC_QUOTIENT`, not evidence that the rational
quotient does not exist.

## Reuse from the typical implementation

Reusable interfaces/oracles in `genus2isogenies.py` and `verify.m` are:

- queue/done graph closure and deterministic vertex sorting;
- modular-Igusa reconstruction, known-model reuse, and exact invariant
  round-trip checks;
- quadratic-twist pinning by Frobenius polynomials;
- reduced model production;
- equality of Frobenius traces, Richelot closure, and period-homomorphism
  degree checks.

They are only partial oracles here: trace equality verifies the underlying
isogeny class, not the principal polarization or RM-marked kernel.  The
existing vertex type is just a curve and loses both `lambda` and `iota`; the RM
BFS needs an internal marked-PPAS vertex.  The existing integer `ell`/1-step/
2-step modular-equation wrappers do not represent an ideal request or a
polarization-only unit transition and must not be reused as though they did.

Phase-1 reducibility tests supply candidate rational primes for unexpected
torsion subrepresentations.  Hurwitz--Maass ideals and unit classes come from
maximal-order class/narrow-sign arithmetic and do not require appearing in
that prime list.  The prime support of the chosen canonical ideal norm should
be tagged `HURWITZ_MAASS_STRUCTURE`, not merged into `reducible_ell`.
Nonmaximal conductor primes belong to the excluded ascent/descent interface.

## Theorem boundary and smallest next slice

Exact and implementable now: field/order validation, class and unit
enumeration, direct (HM) tests, canonical serialization, and all fixture
counts above.

Conditional pending the mathematical/primary-source audits:

- `A[p]` is maximal isotropic in `A[alpha]` over the base field;
- the scaling rule gives exactly the intended isogeny equivalence;
- `U+/U^2` classifies all base-field principal polarizations under the stated
  hypotheses;
- arithmetic pair classes give the complete Hurwitz--Maass action relevant to
  the PPAS problem.

Not part of this slice: nonmaximal ascent/descent, `Q`, `L1`, `L2`, split
cyclic/torsion/inert factorization, termination/completeness of the full note
algorithm, or analytic quotient construction.  In particular, the proof of
the note's full algorithm and its computation of `Q,L1,L2` remain gaps and are
not made true by this prototype.

The recommended next implementation is exactly the CLI above plus byte-stable
JSON and the seven fixtures.  It is useful even if the geometric theorem audit
blocks: in that case it remains an independently checkable arithmetic
reproducer for the disputed class-group condition rather than a speculative
RM graph engine.
