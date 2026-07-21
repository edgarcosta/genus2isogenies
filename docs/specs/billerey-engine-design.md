# Billerey-over-K elliptic engine: design spec

Date: 2026-07-20.

## Context and goal

For an abelian surface A/Q with nontrivial geometric
endomorphisms, compute the primes (or prime ideals) where the mod-ell Galois
representation has nontrivial submodules. Five of the eight endomorphism cases
(QxQ, M2(Q), QxCM, CMxCM, M2(CM)) reduce to questions about elliptic curves over a
number field K: the isogeny primes of each elliptic factor, and the congruence primes
between two factors. This slice builds that reusable elliptic-curve engine in Magma.
No genus-2 wiring yet.

The engine's outputs are statements about the factors as Gal(Qbar/K)-modules; the
surface-level GQ conclusions are the (deferred) case handlers' job. Safe direction:
a GQ-submodule of A[ell] restricts to a GK-submodule. Conditional lemma (quantifiers
matter): if E1[ell] and E2[ell] are both irreducible GK-modules and not isomorphic,
every GK-submodule of E1[ell] + E2[ell] is a coordinate summand; hence a
non-coordinate submodule forces factor reducibility or E1[ell] isomorphic to E2[ell]
(equivalently, in the irreducible case, isomorphic semisimplifications). This yields
candidate supersets; it does not by itself produce an isotropic kernel or a surface
isogeny.

## Public interface

New attachable package `magma/` (with `magma/spec`), file `magma/IsogenyPrimes.m`.
Any implementation MUST match these signatures and semantics exactly.

```
declare verbose IsogenyPrimes, 2;

IsogenyPrimesInfo := recformat<
    Source : MonStgElt,     // always "IsogenyPrimes"
    Kind : MonStgElt,       // "Finite" | "CMFamily"
    Exact : BoolElt,        // certification flag; see Denotational semantics
    IsCM : BoolElt,
    CMOrderDiscriminant : RngIntElt,       // assigned iff IsCM (f^2 * D_F)
    CMFundamentalDiscriminant : RngIntElt, // assigned iff IsCM (D_F)
    CMConductor : RngIntElt,               // assigned iff IsCM (f)
    CMInBaseField : BoolElt,               // assigned iff IsCM
    Stabilized : BoolElt,
    BoundsUsed : SeqEnum >;

CongruencePrimesInfo := recformat<
    Source : MonStgElt,     // always "CongruencePrimes"
    Kind : MonStgElt,       // "Finite" | "AllPrimes" | "Undecided"
    Exact : BoolElt,
    Stabilized : BoolElt,
    BoundsUsed : SeqEnum,
    CertificationMethod : MonStgElt >; // "IsIsogenous" | "ExplicitIsogeny" | "Supplied" | "None"

intrinsic IsogenyPrimes(E::CrvEll :
    AuxBound := 200, MaxAuxBound := 1600,
    FilterBound := 1000) -> SeqEnum, Rec

intrinsic CongruencePrimes(E1::CrvEll, E2::CrvEll :
    NormBound := 1000, MaxNormBound := 8000,
    KnownIsogenous := false,
    CertificationPrimeBound := 100, CertificationDepth := 3) -> SeqEnum, Rec

intrinsic MayBeReducible(ell::RngIntElt, primes::SeqEnum, info::Rec) -> BoolElt
intrinsic MayBeCongruent(ell::RngIntElt, primes::SeqEnum, info::Rec) -> BoolElt
```

Parameter preconditions (`require`d): 0 < AuxBound <= MaxAuxBound,
0 < NormBound <= MaxNormBound, FilterBound > 0, CertificationPrimeBound >= 2,
CertificationDepth >= 0. Base ring must be `FldRat` or `FldNum`. The `Source`
field is set by the producing intrinsic; each accessor `require`s the matching
`Source` (so an isogeny record cannot be fed to `MayBeCongruent` or vice versa)
and `require`s `IsPrime(ell)`. Fields marked "assigned iff IsCM" need `assigned`
guards; the records are the public contract and changes to them are interface
changes.

Dispatch is on ABSOLUTE DEGREE, not Magma type: degree-one fields (including
`RationalsAsNumberField()`) take the exact rational branch. Billerey's theorems
carry extra restrictions at d = 1, so the Billerey branch must never see a
degree-one field. No `printf`; use `vprint IsogenyPrimes`.

Model normalization (first step of both intrinsics): replace each input curve by a
global integral (minimal where Magma supports it) model; every "good/bad" decision
is based on the conductor support of that model, never on the discriminant of the
submitted equation.

### Pinned isogeny primitive (shared helper)

Verified on Magma 2.29-8: `IsogeniesPrimeDegree` and `Isogenies` do not exist;
`IsogenousCurves` and `IsIsogenous` error on CrvEll over FldNum. The one working
primitive over both Q and number fields is `IsogenyFromKernel(E, kernelpoly)`.
Normative helper used by branch 1 and the congruence BFS:

    HasPrimeIsogeny(E, p):
      p = 2: kernel candidates are the linear factors of the 2-division
        polynomial (rational 2-torsion x-coordinates); each linear factor is a
        candidate kernel polynomial.
      p odd: factor the p-division polynomial over K; candidate kernel
        polynomials are products of Galois-stable factor subsets of total degree
        (p-1)/2.
      A candidate is VALIDATED iff `IsogenyFromKernel(E, candidate)` succeeds
      (it errors "Does not appear to be a kernel" otherwise; it is the
      validator). HasPrimeIsogeny is true iff some candidate validates; callers
      may request the FULL LIST of validated codomains from the same candidate
      loop. The certification BFS requires all of them per prime degree; a
      single-witness expansion is incomplete and could return Undecided where
      AllPrimes was reachable.

### Denotational semantics

Let R(E) := { ell prime : E[ell] is reducible as a Gal(Qbar/K)-module }. Let L be
the returned sequence and F(D_F, f) := { p : p split or ramified in the CM field,
p not dividing f }.

| Kind        | guarantee                                                        |
|-------------|------------------------------------------------------------------|
| Finite      | R(E) is a subset of L                                            |
| CMFamily    | F(D_F,f) is a subset of R(E) is a subset of F(D_F,f) union L     |
| AllPrimes   | (congruence) target set = all primes; L is empty by convention   |
| Undecided   | no nontrivial upper bound certified; the only safe upper bound   |
|             | is ALL primes; L is empty; callers MUST short-circuit on Kind    |

`Exact` is a CERTIFICATION flag, not a factual biconditional: if `Exact` is true,
the guarantee upgrades to equality (L = R(E) for "Finite"; the all-primes claim is
proven for "AllPrimes"); if `Exact` is false, NO equality claim is made in either
direction (a candidate list may coincidentally equal R(E)). Assignments: true for
the degree-one branch and for certified "AllPrimes"; false for branch 2, both CM
branches, "Undecided", and uncertified finite congruence output.

For `Kind = "Finite"` the single sequence L serves both phase-2 seeding and the
correctness superset (spanning primes lie in R(E) which lies in L). For
`"CMFamily"` the isogeny class is infinite and seeding the family needs F as well;
the deferred CM-in-K surface handlers consume the record, not L alone.

Accessor formulas (normative, executable):

    MayBeReducible(ell, L, info):
      Kind "Finite":   ell in L
      Kind "CMFamily": (ell does not divide f and KroneckerSymbol(D_F, ell) ne -1)
                       or ell in L
    MayBeCongruent(ell, L, info):
      Kind "Finite":    ell in L
      Kind "AllPrimes": true
      Kind "Undecided": true

### IsogenyPrimes semantics, by branch

The info record always reports geometric CM status via Magma's
`HasComplexMultiplication` (verified: returns geometric CM plus the ORDER
discriminant, over FldRat and FldNum, on Magma 2.29-8). Derived fields:
`CMFundamentalDiscriminant`, `CMConductor`. Normative definition: `CMInBaseField`
is true iff all geometric endomorphisms are defined over K, equivalently
(characteristic zero) iff the CM field embeds in K, tested exactly as
`IsSquare(K ! D_F)`. `HasComplexMultiplication` is analytic; if it fails or is
inconclusive, error loudly (never guess).

1. Absolute degree 1 (checked first, CM or not): `Kind := "Finite"`,
   `Exact := true`, `Stabilized := true`, `BoundsUsed := []`. Pinned algorithm:
   L := { p in {2,3,5,7,11,13,17,19,37,43,67,163} (Mazur's list) :
   HasPrimeIsogeny(E, p) }. Bound parameters ignored. (Empirical note:
   over Q a curve's prime-isogeny set equals its class's, verified to conductor
   400; no class-level pre-pass exists in this branch.)
2. Degree >= 2, no geometric CM: certified superset of R(E); `Kind := "Finite"`,
   `Exact := false`. Every excluded prime is provably irreducible; included primes
   are candidates (phase 2 confirms by constructing isogenies). Local-global
   failures are real (the Sage-documented curve with j = 2268945/128: 7 passes
   every local test over Q and over Q(i) with no global 7-isogeny; verified), so
   `Exact` stays false even after filtering.
3. Degree >= 2, geometric CM: `Kind := "CMFamily"` if `CMInBaseField` else
   `"Finite"`; `Exact := false`; `Stabilized := true`; `BoundsUsed := []`.
   Normative construction: let C_CM := the prime support of the output of Sage
   10.8 `sage.schemes.elliptic_curves.isogeny_class.isogeny_degrees_cm(E)`
   (ported), and

       L := C_CM union { p : p divides f * Disc(K) * Norm(cond E) }.

   Contract lemmas this construction must satisfy (proof sketch normative, full
   write-up deferred to the paper):
   - CM in K (CMInBaseField true): for ell not dividing f * Disc(K) * N(cond E),
     the mod-ell image lies in the Cartan attached to the order O: ell split in O
     gives two stable eigenlines (reducible); ell ramified gives the stable
     kernel line of the nilpotent part (reducible); ell inert gives a nonsplit
     Cartan acting irreducibly. Hence R(E) minus F(D_F, f) is contained in the
     finite fuzz, which is contained in L, giving the CMFamily row of the
     denotation table; the lower containment F(D_F, f) subset R(E) holds by the
     split/ramified eigenline construction.
   - CM not in K (CMInBaseField false): the image lies in the normalizer of the
     Cartan; R(E) is finite and contained in C_CM union the fuzz, i.e. in L
     (Finite row). This is the documented guarantee of `isogeny_degrees_cm`
     (class-group enumeration over the orders between O and the maximal order).

### Non-CM algorithm over K (branch 2)

Sources of truth: Sage 10.8
`sage/schemes/elliptic_curves/gal_reps_number_field.py` (`Billerey_P_l`,
`Billerey_B_l`, `Billerey_R_q`, `Billerey_B_bound`, `Billerey_R_bound`,
`Frobenius_filter`, `reducible_primes_Billerey`) and Billerey, arXiv:0908.1084,
Theorems 2.4 (B_ell) and 2.8 (R_q).

Charpoly helper (the EulerFactor pin was wrong at inert primes,
where `EulerFactor(E, P)` returns a degree-2f polynomial in the rational-prime
variable; verified, and the fixed helper matches Sage's `frobenius_polynomial`):

    Ered := Reduction(E, P);            // good ideals only
    a := TraceOfFrobenius(Ered);
    charpoly := x^2 - a*x + Norm(P);

`Reverse(Coefficients(EulerFactor(E, P)))` may be used only as a cross-check at
residue-degree-1 primes.

Pinned behavior (determinism: output is a pure function of curve and parameters,
scoped to the Magma version recorded in the generated test header; regenerate
expectations on version change, since `Decomposition` order and `IsPrincipal`
generators may vary):

- Auxiliary-prime admissibility: auxiliary rational primes r run over primes
  r >= 5 with r NOT dividing 6 * Disc(K) * Norm(cond E), in increasing order,
  in BOTH phases. Note: Sage's reference skips by the norm of the
  MODEL discriminant, a superset of this; the spec's conductor-based rule is a
  deliberate valid tightening and is the second pinned deviation from Sage.
  B_r is Billerey's B_ell
  exactly (Theorem 2.4): the composed product per eq (9) runs over ALL primes q
  above r with the 12*e_q exponents and residue-degree brackets.
- Zero-gcd sentinel: while an accumulator is 0, its candidate prime set is a
  formal TOP element that never satisfies "unchanged, stop" and is never passed
  to `PrimeFactors` (verified: `PrimeFactors(0)` errors in Magma and Sage).
  `PrimeFactors` only on nonzero final accumulated gcds.
- B-phase with auto-escalation: accumulate B := GCD(B, B_r) over admissible r up
  to the current bound, starting at AuxBound. After exhausting a bound, compare
  the candidate prime set with the previous bound's: unchanged across a doubling
  (including a comparison landing exactly at the cap): stop with
  `Stabilized := true`; else double, capped at MaxAuxBound; cap reached while
  still shrinking: `Stabilized := false`. When AuxBound = MaxAuxBound (a single
  evaluation, no doubling possible), `Stabilized := (B eq 1)`. Early exit any
  time B reaches 1 (absorbing; counts as stabilized). The plateau rule is a
  reproducibility/tightness heuristic, not a convergence certificate; safe
  because every intermediate set is already a certified superset.
- If the B-phase ends with B still 0, run the R-phase over prime ideals q above
  admissible rational primes, in increasing norm (within a rational prime,
  Magma's `Decomposition` order), same escalation rule on the norm bound.
  Reality check: Billerey's Example 5.8 curve over
  K = Q(sqrt -3, sqrt -7) (conductor 2*O_K, non-CM) has B_l = 0 for every
  admissible l, so this phase genuinely runs. Pinned formula: the FULL Sage
  `Billerey_R_q`, product over k = 0..d/2 INCLUDING the k = 0 term (which
  contributes P(1)).
- Lemma (normative): Billerey's Theorem 2.8 proof uses only that q^{h_K} is
  principal, through a generator gamma whose minimal polynomial enters the
  resultants; replacing h_K by any h with q^h principal, in particular
  h = order of [q] in Cl(K), preserves the conclusion verbatim with gamma a
  generator of q^h. No new exceptional primes arise. This is the first pinned
  deviation from Sage. Generator = whatever `IsPrincipal` returns; R_q integers
  are generator-dependent up to units, so comparison is on FINAL outputs only.
- Gate: for l INERT in K and q = (l), R_q = B_l; asserted for inert primes
  only. For split principal q the equality FAILS (B_l can vanish while R_q does
  not; verified over Q(sqrt 29), and Sage's own docstring shows B_l = 0 at split
  and nonzero at inert primes for that curve).
- If the cap is reached with both phases' gcds still 0: error (Billerey
  guarantees nonvanishing R_q for some q on non-CM curves; persistent 0 means CM
  was missed or the cap is absurdly small).
- Candidate set := prime factors of the final gcd, union {2, 3}, union prime
  factors of Disc(K), union residue characteristics of the conductor support of
  E. A deliberate conservative ENLARGEMENT of the theorems' excluded set;
  tightness soft-compared against Sage.
- Frobenius filter, applied once after escalation: discard candidates ell having
  a good ideal q, Norm(q) <= FilterBound, whose Frobenius charpoly is irreducible
  mod ell with constant term nonzero mod ell.
- Return sorted survivors; `BoundsUsed := [ final_B_bound, final_R_bound_or_0,
  FilterBound ]`.

### CongruencePrimes semantics

Target: all primes ell such that E1[ell] and E2[ell] have isomorphic
semisimplifications over K. Output per the denotation table; `Exact` true only
for certified `"AllPrimes"`.

Flow:

0. If `KnownIsogenous` (consumer-supplied relation, e.g. the M2 handlers):
   return `"AllPrimes"`, `Exact := true`,
   `CertificationMethod := "Supplied"`, `Stabilized := true`,
   `BoundsUsed := [0, 0, 0]`. Required because Magma 2.29-8's `IsIsogenous`
   errors on CrvEll over FldNum (verified), so over number fields in-engine
   certification can only confirm by explicit search.
1. Degree 1: call `IsIsogenous` FIRST. Isogenous: `"AllPrimes"`,
   `Exact := true`, `CertificationMethod := "IsIsogenous"`,
   `Stabilized := true`, `BoundsUsed := [0, 0, 0]`. Not isogenous: proceed,
   remembering the answer is finite (but see step 4: a finite cap may still
   fail to separate). Rule for every path: `Stabilized := true` on every
   certified (`Exact`) return; `false` on `"Undecided"`; on `"Finite"` it
   reports the escalation plateau rule.
2. Gcd loop: over prime ideals q of K good for BOTH curves (conductor support
   after model normalization), in increasing norm, accumulate
   G := GCD over q of ( p_q * (a_q(E1) - a_q(E2)) ), p_q the residue
   characteristic. The p_q factor keeps ell = p_q admissible at its own q; with
   bad places skipped this covers ell = 2, quadratic twists, and bad primes.
   Zero-sentinel and escalation exactly as in the B-phase (NormBound to
   MaxNormBound; the single-evaluation Stabilized rule applies).
3. If G != 0: return sorted prime factors of G, `Kind := "Finite"`,
   `Exact := false`, `CertificationMethod := "None"`,
   `BoundsUsed := [final_norm_bound, 0, 0]`.
4. If G = 0 at MaxNormBound (traces agreed everywhere):
   - Degree 1 with `IsIsogenous` false (this state IS reachable;
     a finite cap does not guarantee a separating prime below it): return
     `"Undecided"`, `Exact := false`, `Stabilized := false`,
     `CertificationMethod := "None"`, `BoundsUsed := [MaxNormBound, 0, 0]`.
   - Degree >= 2: bounded explicit-isogeny certification: BFS from E1 through
     prime-degree isogenies built by `HasPrimeIsogeny` (validated kernels) of
     degree <= CertificationPrimeBound, depth <= CertificationDepth. Node match
     test (twist soundness): a node matches iff
     `IsIsomorphic(node, E2)` over K. Equal j-invariants are ONLY a prefilter
     before calling `IsIsomorphic`; bare j-equality must never certify (a
     nontrivial quadratic twist has the same j but twists the representation).
     Deduplicate BFS nodes by K-isomorphism within equal-j buckets. Found:
     `"AllPrimes"`, `Exact := true`, `CertificationMethod := "ExplicitIsogeny"`,
     `Stabilized := true`,
     `BoundsUsed := [MaxNormBound, CertificationPrimeBound, CertificationDepth]`.
     Not found: `"Undecided"`, `Exact := false`, `Stabilized := false`, same
     BoundsUsed (search can only confirm, never refute).

No further filtering of the finite output in this slice.

## Dependencies

Magma >= 2.28 with CHIMP attached for the charpoly star-product intrinsics. The
package routes CHIMP intrinsics through dynamic by-name lookup: `BillereyBl`,
`BillereyRq`, and the IsogenyPrimes Billerey branch carry functional
attach-CHIMP `require`s, while the degree-one (over-Q) paths run without CHIMP.
SageMath (10.8 here) only for generating test data. NOT used anywhere: PARI
`ellisomat` (broken over number fields in the installed PARI), Magma
`IsIsogenous`/`IsogenousCurves` over FldNum, `IsogeniesPrimeDegree`, `Isogenies`
(absent on 2.29-8).

## Error handling, style, implementation notes

- `require` on preconditions (base ring type, common base field after
  normalization, parameter ordering, accessor Source/prime checks).
- `vprint IsogenyPrimes, 1` (progress) and `, 2` (per-auxiliary detail: bound,
  gcd support, stopping reason).
- Zero-gcd sentinel as pinned; `PrimeFactors` only on nonzero final gcds.
- Degree/bit-size sanity caps with loud diagnostics; cache Frobenius data and
  ideal decompositions per curve.
- Records and accessors are the public contract; document fields, `assigned`
  discipline, and the denotation table at the declaration site.

## Risks and known subtleties

- Potential CM only over an extension: branch selection tests geometric CM;
  `HasComplexMultiplication` is analytic and must fail loudly if inconclusive.
- The plateau rule can stop before the minimal superset; `Stabilized` and
  `BoundsUsed` expose that.
- Version-scoped determinism: `Decomposition` order and `IsPrincipal` generators
  can differ across Magma versions; expectations record the version and must be
  regenerated on upgrade.
- Non-CM output is a superset even after filtering; consumers use the accessors,
  never raw membership.
- CertificationDepth/PrimeBound bound the isogeny search; a true isogeny of
  larger degree yields Undecided (safe; the M2 handlers should pass
  KnownIsogenous instead of relying on the search).
- Kernel-polynomial enumeration in `HasPrimeIsogeny` is exponential in the
  number of division-polynomial factors in principle; harmless at Mazur/BFS
  sizes, capped loudly if a pathological case appears.

