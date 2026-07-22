# Wave 3 adversarial review

Date: 2026-07-22 UTC

## Verdict

**No blocker was found in the frozen maximal-order arithmetic slice or in the
paper proposition.**  The implementation distinguishes the two conditions

\[
  \mathfrak a^2\text{ is principal}
  \quad\text{and}\quad
  \mathfrak a^2=(\alpha)\text{ for some }\alpha\gg0,
\]

enumerates the latter condition directly, and does not turn arithmetic
Hurwitz--Maass elements into claims about computed quotient PPASs.

This verdict is deliberately split in two:

- **Arithmetic validation:** exact ideal, unit, sign, norm, and pair-equivalence
  attacks passed for the named fields below.  The public suite passed all 13
  tests.  The adversarial reconstruction uses a different control flow but
  shares Sage's class-group, ideal, and unit backends, so it is corroborating
  exact evidence rather than a fully independent computer-algebra oracle.
- **Geometric nonvalidation:** no abelian surface, polarization, RM action,
  quotient, descent datum, or polarized-isomorphism certificate is input or
  computed.  No result here validates target existence or uniqueness, the
  analytic quotient layer, or the full RM graph.

Disposition words below have their literal meanings: **PASS** is a successful
attack in the frozen scope; **ACCEPTED boundary** is an intentional limit that
must survive into a production API; **GAP/OUT_OF_SCOPE** is not established by
this slice; and **BLOCKER** would invalidate a promised scoped result.  There
are no findings with disposition **BLOCKER**.

## Arithmetic attacks and dispositions

| ID | Attack and exact evidence | Disposition |
|---|---|---|
| AR-01 | **Duplicate classes in `U^+/U^2`.**  For `D=12`, the output is `1` and `2+w`.  The latter is totally positive and has odd free-unit coordinate, so their ratio is not a unit square.  For `D=60`, two unit classes times two eligible ideal classes give four records; an exact pairwise equivalence test found no equivalent pair. | **PASS** |
| AR-02 | **Missing classes in `ker(Sq)`.**  Directly enumerating every ordinary class up to the real-quadratic Minkowski bound and testing its square for a positive generator gave exactly the output coordinates: `D=40`: `(0),(1)`; `D=145`: `(0),(2)` inside ordinary `C4`; `D=205`: `(0)`; `D=229`: `(0)`.  No eligible class was omitted in the seven acceptance fields. | **PASS** (same Sage arithmetic backend; not a fully independent CAS result) |
| AR-03 | **Assigning an `alpha` to a nonprincipal ideal.**  In `D=40`, the selected norm-two ideal has nonzero ordinary class coordinate `(1)`, so the ideal itself is nonprincipal.  Its square, not the ideal, is `(2)`.  The record correctly says `a^2=(alpha)` with `alpha=2`; it never says `a=(alpha)`. | **PASS** |
| AR-04 | **Confusing ordinary-principal with narrow-principal.**  In `D=205`, the nontrivial ordinary class has canonical representative of norm `3`, row HNF `[[1,2],[0,3]]`, and square generator `w-7`.  The two signs of `w-7` are opposite.  All four unit-parity multiples still have mixed signs because the unit-signature image is only `{(-1,-1),(1,1)}`.  The implementation therefore excludes this ordinary two-class even though its square is principal. | **PASS** |
| AR-05 | **Using the image rather than the kernel of squaring.**  In `D=229`, both ordinary and narrow groups are `C3`.  Squaring is surjective on `C3` but has trivial kernel.  The output has one ideal class and one record, so it does not include the two false image-of-squaring candidates. | **PASS** |
| AR-06 | **Missing an eligible class when the image is too small.**  In `D=40` (the field `Q(sqrt(10))`), the ordinary `C2` class with norm-two representative is retained even though `Sq:Cl(F)->Cl^+(F)` has trivial image. | **PASS** |
| AR-07 | **Nonmaximal order and order-discriminant confusion.**  `D=5,f=2` exits `2` with `UNSUPPORTED_NONMAXIMAL_ORDER`, order discriminant `20`, and no HM data.  `D=20,f=1` exits `2` with `INVALID_FIELD_DISCRIMINANT`; it is not silently reinterpreted as the conductor-two order in the field of discriminant `5`.  Conductors `0` and `-1` return `INVALID_ORDER_CONDUCTOR`. | **PASS** |
| AR-08 | **Exact record reconstruction.**  For `D=5,12,40,60,145,205,229`, reconstructing every output ideal from its row HNF and every `alpha` from its `[1,w]` coefficients verified `a^2=(alpha)` and positivity at both real embeddings.  An all-pairs test using principal ideal quotients and exact unit-group coordinates found no two equivalent output pairs. | **PASS** (corroborating Sage computation) |
| AR-09 | **Canonicalization and identity drift.**  The successful `D=60` process output was byte-identical on repetition, and the public test checks one JSON line, empty stderr, and a final newline.  A stress sweep through all 153 positive fundamental discriminants `D<=500` found no failed run check, wrong genus-theory count, or case where `unit-0` was not `[1,0]`.  Ideal selection is by the promised `(norm,row-HNF)` key. | **PASS** for the fixed Sage toolchain; **ACCEPTED boundary** across Sage versions because diagnostic class coordinates depend on Sage's chosen class-group generators and are not claimed as a cross-version interchange format. |
| AR-10 | **Dependence on the field involution/chosen RM embedding.**  Output coefficients use the fixed canonical generator `w`; no quotient by `w |-> Tr(w)-w` is added.  For the nontrivial eligible class of `D=145`, the chosen ideal has HNF `[[1,3],[0,4]]`, while its conjugate has HNF `[[4,0],[0,1]]`.  Their ideal quotient is principal and the induced positive unit has unit coordinates `[0,-2]`, a square, so this example returns to the same arithmetic pair through the declared pair equivalence, not through an illicit field-automorphism quotient.  A future surface adapter must still interpret the coefficients through the supplied `iota`; conjugating `iota` is different marked input. | **PASS / ACCEPTED boundary** |
| AR-11 | **Validated error/status behavior.**  No command, an unknown command, a noninteger discriminant, and an extra flag each exit `2`, emit one JSON `INVALID_ARGUMENTS` document, and put the human diagnostic on stderr.  Mathematical input rejections contain no HM records.  Generic proof/internal exceptions are statically routed to exit `3` statuses rather than nonexistence. | **PASS** for exercised success and rejection paths; exit-`3` fault injection is **NOT_INDEPENDENTLY_VALIDATED**. |
| AR-12 | **Help output versus the one-JSON wording.**  `--help` and `enumerate --help` are argparse metacommands: they exit `0` with 12 and 7 lines of text, respectively, rather than JSON.  They are outside the frozen `enumerate --field-discriminant D --order-conductor f` result grammar, but the unqualified sentence “Stdout is one byte-stable JSON document” could be read more broadly. | **GAP/OUT_OF_SCOPE**, nonblocking.  A production contract should explicitly exempt CLI help or make help a JSON status. |
| AR-13 | **Status-specific rejection metadata.**  The nonmaximal rejection correctly has the decisive status and no HM records, but its shared `theorem_boundary` says `arithmetic_result=EXACT_PROOF_ENABLED` even though class/unit arithmetic was intentionally not run. | **GAP/OUT_OF_SCOPE**, nonblocking semantic overbreadth.  Production result variants should have status-specific payloads. |
| AR-14 | **Unbounded normalization.**  Unit-square normalization has a defensive 10,000-step cap.  No input bound is stated, and no proof of that operational cap was found.  Exceeding it would, however, return an explicit exit-`3` internal-inconsistency state rather than a false arithmetic answer. | **ACCEPTED boundary** for a throwaway prototype; the production implementation should use a directly bounded/logarithmic reduction or document a resource limit. |

The sweep in AR-09 used the independent classical cardinality formula
`2^(omega(D)-1)` as its expected value.  That formula is not called by
`prototype.sage`; nevertheless, the field and class computations in the sweep
still ran in Sage, so the result is not presented as an independent CAS check.

## Geometric, semantic, and paper attacks

| ID | Attack and evidence | Disposition |
|---|---|---|
| GEO-01 | **Confusing the underlying abelian variety with a PPAS.**  The frozen geometric input is `(A,lambda,iota)`, while the prototype takes only `(D,f)`.  Success is named `OK_ARITHMETIC_ONLY`; top-level `surface_scope` is `NOT_EVALUATED`; and the theorem boundary says both quotient PPAS and target deduplication are not computed. | **PASS** |
| GEO-02 | **Treating an arithmetic HM element as a unique target PPAS.**  Every record labels its geometric degree/weight interpretation `THEOREM_CONDITIONAL`.  The report explicitly warns about stabilizers and requires exact polarized `k`-isomorphism deduplication.  Thus record distinctness is only distinctness in the arithmetic pair quotient. | **PASS** |
| GEO-03 | **Field-of-definition slippage.**  The scoped statement uses `iota:O_F -> End_k(A)`, a reference principal polarization over `k`, and polarized `k`-isomorphism.  It does not infer that geometric endomorphisms or geometric isomorphisms descend to `k`. | **PASS** |
| GEO-04 | **Confusing `End_k` with `End_kbar`.**  The JSON calls larger geometric endomorphisms out of scope.  The report and paper proposition both state that `End_kbar(A)` may be larger.  The changed Rosati formula is asserted on `End_kbar^0(A)`, but only its restriction to the commutative `End_k^0(A)=iota(F)` is asserted to remain the identity. | **PASS** |
| GEO-05 | **Forgetting the chosen embedding.**  The paper fixes an actual isomorphism into `End_k(A)`, not merely an abstract equality of rings.  The future analytic return contract requires the transported embedding.  The arithmetic prototype fixes a field presentation but correctly cannot certify a surface action that it never receives. | **PASS / ACCEPTED boundary** |
| GEO-06 | **Paper proposition versus cited hypotheses.**  González--Guàrdia--Rotger, Theorem 2.10, assumes an abelian variety of `GL_2`-type over `k`, sets `R=End_k(A)`, and assumes a reference principal polarization over `k`.  In the real case it gives positive units modulo unit squares.  The insertion assumes the stronger concrete specialization `dim(A)=2` and `End_k(A)=iota(O_F)` for maximal real quadratic `O_F`, so all cited hypotheses are present.  Geometric simplicity is an additional restriction, not an extrapolation. | **PASS** ([primary article, DOI 10.4064/aa116-3-3](https://doi.org/10.4064/aa116-3-3)) |
| GEO-07 | **Uniqueness and automorphism action in the paper.**  Relative to the fixed principal polarization, the Neron--Severi/symmetric-endomorphism correspondence gives a unique positive unit.  Since `Aut_k(A)=iota(O_F^x)` and Rosati is the identity on the real field, pullback acts by `epsilon |-> u^2 epsilon`.  The proposition says “on the fixed underlying surface” before passing to polarized isomorphism, so it does not confuse uniqueness of a representative with uniqueness of a polarized class. | **PASS** |
| GEO-08 | **Changing the reference polarization.**  Direct substitution in the Rosati definition gives `x^(dagger_eta)=iota(eta)^-1 x^dagger iota(eta)`.  Commutativity makes the restriction to `iota(F)` unchanged; no analogous claim is made for extra geometric endomorphisms. | **PASS** |
| GEO-09 | **Unsupported geometric cases.**  Nonmaximal orders, RM not defined over `k`, larger rational endomorphism rings, products, QM, matrix algebras, and geometric polarization classification are explicitly excluded.  The paper diff contains only the fixed-surface proposition and no Hurwitz--Maass completeness or graph claim. | **PASS / ACCEPTED boundary** |
| GEO-10 | **Analytic quotient and target equality.**  No exact quotient, rational or `k`-model recognition, twist/descent, transported action, kernel certificate, or polarized-isomorphism key was exercised in this wave.  Approximate modular invariants would not supply these certificates. | **GAP/OUT_OF_SCOPE**; **NOT GEOMETRICALLY VALIDATED**. |

The paper build result is preserved in `PAPER_LATEXMK.md`; this adversarial pass
inspected the diff and hypotheses but did not treat that existing build
transcript as an independent mathematical validation.

## Exact counterexamples retained

These examples are important because replacing the direct predicate by a
nearby class-group slogan changes the answer.

1. **`D=40`: kernel versus image.**  The nontrivial ordinary `C2` class has
   norm-two representative `a`, with `a^2=(2)` and `2>>0`.  It must be
   included even though `Sq:Cl(F)->Cl^+(F)` has trivial image.
2. **`D=205`: principal versus narrow-principal.**  The nontrivial ordinary
   `C2` class has `a^2=(w-7)`, but `w-7` and every unit-parity multiple have
   mixed signs.  It must be excluded.
3. **`D=229`: kernel versus surjective image.**  Squaring on `C3` is an
   automorphism, but only the identity lies in its kernel.  An image-based
   enumeration would add two false classes.
4. **`D=5,f=2` versus `D=20,f=1`: order versus field discriminant.**  The first
   is a valid field with an unsupported nonmaximal order; the second is not a
   fundamental field discriminant.  The two inputs require different statuses.

## Commands and observed summaries

The complete public-boundary suite was rerun as

```sh
PYTHONHASHSEED=0 sage -python docs/superpowers/rm-goal/test_prototype.py -v
```

and ended with:

```text
Ran 13 tests in 22.870s
OK
```

The exact adversarial reconstruction of eligible classes and pairwise
equivalence produced:

```text
D=5 eligible=[()] output=[()] records=1 exact=True equivalent-pair-duplicates=[]
D=12 eligible=[()] output=[()] records=2 exact=True equivalent-pair-duplicates=[]
D=40 eligible=[(0,), (1,)] output=[(0,), (1,)] records=2 exact=True equivalent-pair-duplicates=[]
D=60 eligible=[(0,), (1,)] output=[(0,), (1,)] records=4 exact=True equivalent-pair-duplicates=[]
D=145 eligible=[(0,), (2,)] output=[(0,), (2,)] records=2 exact=True equivalent-pair-duplicates=[]
D=205 eligible=[(0,)] output=[(0,)] records=2 exact=True equivalent-pair-duplicates=[]
D=229 eligible=[(0,)] output=[(0,)] records=1 exact=True equivalent-pair-duplicates=[]
```

The additional bounded stress command loaded `prototype.sage`, ran every
positive fundamental discriminant through `500`, checked all run checks and
the independent genus-theory cardinality, and printed:

```text
fundamental_fields 153 max_records (4, 492) failures []
```

These two auxiliary checks are exactly reproducible from Sage objects, but
because they share Sage arithmetic with the prototype they are explicitly
labelled **NOT FULLY INDEPENDENT**.

## Unresolved broader gaps

The adversarial review found that `REPORT.md` keeps the following claims out of
the positive slice.  Their disposition remains **GAP/OUT_OF_SCOPE**, not an
implicit success:

- ramified primes need a third local case; defining “split” only by norm `p`
  incorrectly includes them;
- the note's bound `m<=(e-1)/2` excludes the valid `e=2`, pure-torsion kernel;
- the displayed cyclic-factor subgroup uses the wrong divisor exponent and
  fails already at `e=1`;
- finite effective construction of the node-dependent `Q,L_1,L_2` is absent;
- Faltings finiteness does not by itself give executable termination,
  decidable target equality, or completeness; and
- the analytic reconstruction layer does not yet satisfy the exact return
  contract.

None of these gaps affects the truth of the emitted maximal-order arithmetic
pairs.  Each blocks expansion from those requests to the full RM graph until
its own lemma, algorithm, and certificate are supplied.

## Final disposition

The frozen prototype and the scoped paper proposition may survive this review.
The safe handoff is the arithmetic list of pair classes together with explicit
surface/analytic obligations.  It is not a computed isogeny graph.  A
production extraction should preserve the rejection variants as distinct
types, make their payloads status-specific, keep the chosen `iota` in every
geometric call, and require an exact polarized `k`-isomorphism certificate
before deduplicating targets.
