# Frozen Wave 2 slice: maximal-RM arithmetic only

This specification was frozen after the three independent Wave 1 audits and before prototype implementation. It is the only implementation target in this run.

## Scope and nonclaims

The geometric consumer is a marked PPAS `(A, lambda, iota)` over a number field `k`, with `A` geometrically simple and a chosen isomorphism `iota: O_F -> End_k(A)`, where `O_F` is the maximal order of a real quadratic field. `End_k(A)`, `End_kbar(A)`, the chosen embedding, and the Rosati-fixed subring are distinct data. The prototype itself receives no surface and therefore returns arithmetic requests only. It does not compute or deduplicate quotient PPASs, prove the full RM graph algorithm, or support nonmaximal orders, RM not defined over `k`, larger rational endomorphism rings, products, QM, or matrix algebras.

## Public seam

The sole tested interface is:

```sh
PYTHONHASHSEED=0 sage -python docs/superpowers/rm-goal/prototype.sage \
  enumerate --field-discriminant D --order-conductor f
```

`D` must be a positive fundamental real-quadratic discriminant. The canonical field presentation is `O_F = Z[w]`, with `w^2-w+(1-D)/4=0` for odd `D` and `w^2-D/4=0` for even `D`. The supplied order is `Z[f w]`; `f != 1` is rejected, and equality with Sage's maximal order is checked rather than inferred. Proof-enabled class-group, narrow-class-group, principal-ideal, and unit computations are mandatory.

Stdout is one byte-stable JSON document. Exit `0` means `OK_ARITHMETIC_ONLY`, exit `2` is a validated input/scope rejection, and exit `3` is an unknown proof failure or internal inconsistency. Exceptions are never converted into mathematical nonexistence.

## Exact output

Let

`Sq: Cl(F) -> Cl+(F), [a] |-> [a^2]+`.

The program returns:

1. deterministic representatives of `U+/U^2`, where `U = O_F^x` and `U+` is the totally positive subgroup;
2. deterministic ordinary ideal representatives for `ker(Sq)`, selected by least `(norm, row-HNF)` among all integral ideals up to the real-quadratic Minkowski bound `floor(sqrt(D)/2)`; and
3. the duplicate-free arithmetic pair representatives `(a, epsilon*alpha_a)`, where `a^2=(alpha_a)`, `alpha_a` is totally positive, and `epsilon` runs through `U+/U^2`.

This resolves the note's discrepancy by implementing the initial condition `a^2=(alpha), alpha >> 0`. The later phrase about “classes in Cl+(F) that are squares of elements of Cl(F)” is not implemented: it is ill-typed and is neither equivalent to `ker(Sq)` nor a safe proxy for it.

Every record includes the ideal norm/HNF, ordinary class coordinates as diagnostics, exact element coefficients in `[1,w]`, both signs, expected quotient degree/graph weight labelled theorem-conditional, and checks proving `a^2=(alpha)`, `alpha >> 0`, and `Norm(alpha)=Norm(a)^2`. Representatives are normalized modulo square units using exact algebraic-real comparisons. Conjugate ideals are not identified because the chosen RM embedding is part of the future marked input.

The JSON also reports ordinary/narrow class invariants, unit-signature data, the arithmetic-only theorem boundary, and explicit unsupported/unknown state vocabulary. It records, but does not execute, the fixed-surface statement that principal polarizations are `lambda_epsilon=lambda o iota(epsilon)` and are classified up to polarized `k`-isomorphism by `U+/U^2`. Replacing `lambda` by `lambda_eta` changes Rosati to `x |-> eta^-1 x^dagger eta`; its restriction to the commutative `End_k^0(A)=F` stays the identity, while its effect on larger geometric endomorphisms is out of scope.

## Acceptance fixtures and independent checks

Public-boundary tests use literal expectations for discriminants `5, 12, 40, 60, 145, 205, 229`: negative/positive unit norm, unequal ordinary/narrow groups, admissible nontrivial 2-torsion, an ordinary `C4`, an ineligible ordinary 2-class, and the `Cl=C3` image-versus-kernel trap. `D=5,f=2` must return `UNSUPPORTED_NONMAXIMAL_ORDER`; `D=20,f=1` must return `INVALID_FIELD_DISCRIMINANT`. Repeated success output must be byte-identical.

For every fundamental `D` in the brute-force range, the HM-record count is checked independently against the classical Hurwitz--Maaß/genus-theory value `2^(omega(D)-1)` and against the 2-torsion cardinality of the narrow class group. These checks validate cardinality, not geometric target uniqueness.

## Required future analytic result

An RM BFS may consume a pair only when its quotient engine returns either an explicit unsupported/unknown state or a certified target PPAS and polarization, transported `O_F` embedding, exact kernel/degree/pullback identity, rational model with twist/descent data, and a polarized `k`-isomorphism certificate/key. Approximate modular invariants alone are insufficient.
