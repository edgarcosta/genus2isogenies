# Maximal real-quadratic RM: audited algorithm bridge

## Decision

The first safe implementation slice is exact arithmetic for the maximal-order Hurwitz--Maass data and for principal polarizations on a fixed underlying surface. The full RM breadth-first graph in `private/tex/note.tex` is **not** justified: its ramified-prime handling, a local normal form, the cyclic factorization formula, the effective construction of `Q,L_1,L_2`, and executable termination/completeness all have independent gaps.

Positive statements in this report have the following fixed scope:

\[
 (A,\lambda,\iota)/k,\qquad
 \iota:\mathcal O_F\xrightarrow{\sim}\operatorname{End}_k(A),
\]

where `k` is a number field, `A` is geometrically simple, `F` is real quadratic, `O_F` is its maximal order, `iota` is defined over `k`, and `lambda` is a principal polarization over `k`. This does **not** assert `End_kbar(A)=O_F`. A chosen RM embedding, a different conjugate embedding, the rational endomorphism ring, the geometric endomorphism ring, and the Rosati-fixed subspace are distinct data.

Nonmaximal orders, RM not defined over `k`, larger rational endomorphism rings, products, QM, `M_2(Q)`, and CM are preserved as explicit unsupported interfaces. No classification below is extrapolated to them.

## Source audit

The accepted dependencies were checked against primary sources:

- Mumford's *Abelian Varieties*, §20(3), §21 Theorem 1, Application III, and §23 supplies the Néron--Severi/symmetric-endomorphism correspondence, positivity of Rosati, the ample cone, pairing functoriality, descent through isotropic kernels, and the maximal-isotropic order criterion ([primary scan](https://wstein.org/edu/Fall2003/252/references/mumford-abvar/Mumford-Abelian_Varieties.pdf)).
- González--Guàrdia--Rotger, Theorem 2.10, gives the base-field classification of principal polarizations for real `GL_2`-type varieties and its automorphism action ([Acta Arith. 116 (2005), DOI 10.4064/aa116-3-3](https://doi.org/10.4064/aa116-3-3)).
- van der Geer, *Hilbert Modular Surfaces*, §I.4, pp. 12--13, defines `Sq:Cl(F)->Cl+(F)`, chooses representatives of its kernel, and describes the positive-unit kernel of the Hurwitz--Maass extension ([book DOI](https://doi.org/10.1007/978-3-642-61553-5)).
- Hauffe-Waschbüsch--Krieg give the ideal-square description and index `2^(omega(D)-1)` of the classical Hurwitz--Maaß extension in equation (8) ([Research in Number Theory 8 (2022), article 47](https://doi.org/10.1007/s40993-022-00346-5)).
- Faltings proves the relevant finiteness theorem, but not the missing effective neighbor enumeration in the note ([Invent. Math. 73 (1983), DOI 10.1007/BF01388432](https://doi.org/10.1007/BF01388432)).

The local `private/tex/note.tex` was compared once, byte-for-byte, with the official source of [arXiv:2405.19820v1](https://arxiv.org/abs/2405.19820v1). `ARXIV_DIFF.log` records the complete diff. The differences are only abstract wording and one in-preparation bibliography entry; there are no mathematical changes. These are copies of the same work and were not treated as independent validation.

## Exact objects and data layers

### Arithmetic input

The prototype's exact input is a positive fundamental discriminant `D` and an order conductor `f`. It constructs the canonical real quadratic field presentation `O_F=Z[w]` and the supplied order `Z[f w]`. Only `f=1` and equality with Sage's certified maximal order are supported.

Arithmetic output consists of ordinary and narrow class-group diagnostics, totally positive unit square classes, canonical integral ideals, positive generators, and exact equality/sign/norm certificates. It contains no surface, polarization, period matrix, or approximate number.

### Polarized input for the theorem and future adapter

The geometric consumer must additionally supply the marked object `(A,lambda,iota)` above. `lambda` is not decorative: changing it changes Rosati on the whole geometric endomorphism algebra. `iota` determines how ideals act and is retained through every quotient.

### Analytic input

Period matrices, analytic action matrices, interval precision, modular invariants, rational-recognition data, curve/product models, and twist data belong only to a future quotient engine. They are not inferred from the arithmetic prototype.

## Verified finite maximal-order data

Define

\[
 \operatorname{Sq}:\operatorname{Cl}(F)\longrightarrow\operatorname{Cl}^+(F),
 \qquad [\mathfrak a]\longmapsto[\mathfrak a^2]^+.
\]

This map is well-defined because replacing `a` by `x a` changes its square by the narrow-principal ideal `(x^2)`. A Hurwitz--Maass pair is

\[
 (\mathfrak a,\alpha),\qquad \mathfrak a^2=(\alpha),\quad \alpha\gg0,
\]

modulo `(a,alpha)~(x a,x^2 alpha)` for `x in F^x`. Pair classes form an elementary 2-group `H_F` in an exact sequence

\[
 1\longrightarrow
 \mathcal O_F^{\times,+}/(\mathcal O_F^\times)^2
 \longrightarrow \mathcal H_F
 \longrightarrow \ker(\operatorname{Sq})
 \longrightarrow 1.
\]

The proof is direct: projection remembers the ordinary ideal class; a pair in the kernel scales to `(O_F,epsilon)` with `epsilon` a positive unit; two such pairs differ by a unit square. Thus a duplicate-free set of **arithmetic pair classes** is obtained by choosing one canonical ideal `a_c` for each `c in ker(Sq)`, one verified positive generator `alpha_c` of `a_c^2`, and emitting `(a_c,epsilon alpha_c)` for canonical `epsilon` in the positive-unit quotient.

This reconciles the note's two formulations. The initial condition `a^2=(alpha), alpha>>0` is correct. The later phrase “classes in `Cl+(F)` that are squares of elements of `Cl(F)`” is ill-typed and wrong. The prototype therefore tests the initial condition directly. `Q(sqrt(10))` shows that an image-of-squaring rule can miss an eligible class; `Q(sqrt(229))`, where ordinary and narrow class groups are `C3`, shows that it can add ineligible classes.

These representatives have no duplicates as marked arithmetic pair classes. Their action on a special moduli point can have a stabilizer, so target PPASs still require exact polarized `k`-isomorphism deduplication; arithmetic equality alone is not a target-isomorphism oracle.

## Principal polarizations on the fixed surface

Rosati restricts to the identity on `End_k^0(A)=iota(F)`: a positive involution cannot be the nontrivial real-quadratic automorphism. Relative to `lambda`, every `k`-defined principal polarization on the fixed underlying `A` is

\[
 \lambda_\epsilon=\lambda\circ\iota(\epsilon),
 \qquad \epsilon\in\mathcal O_F^{\times,+}.
\]

Since `Aut_k(A)=iota(O_F^x)` and Rosati is the identity there, pullback by `iota(u)` sends `epsilon` to `u^2 epsilon`. Hence the polarized `k`-isomorphism classes are exactly

\[
 \mathcal O_F^{\times,+}/(\mathcal O_F^\times)^2.
\]

For a real quadratic field this set has one element when a negative-norm unit exists and two otherwise. If the reference is changed to `lambda_eta`, then

\[
 x^{\dagger_{\lambda_\eta}}=\iota(\eta)^{-1}x^{\dagger_\lambda}\iota(\eta).
\]

The restriction to `iota(F)` remains the identity; the involution on a larger `End_kbar^0(A)` may change and is out of scope.

## What the analytic quotient engine must return

For RM breadth-first enumeration to continue, a successful quotient call on `(A,lambda,iota;a,alpha)` must return:

1. an exact target PPAS `(B,lambda_B)` over `k` and a target representation/model;
2. the transported embedding `iota_B:O_F -> End_k(B)`, with exact analytic/algebraic action data sufficient to form the next ideal kernel;
3. an exact kernel, ideal, degree, and pullback certificate `phi^vee lambda_B phi = lambda iota(alpha)`;
4. certified rational invariants/model recognition plus twist and descent data;
5. an exact polarized `k`-isomorphism key and certificate for every deduplication, recording whether the RM marking is preserved, conjugated, or forgotten; and
6. precision and trace/hash evidence usable by the existing reconstruction and Frobenius oracles.

Failure to certify any item is an explicit `UNKNOWN_*` or `UNSUPPORTED_*` result. Approximate Igusa invariants are not a negative answer and are not a `k`-isomorphism key.

## Phase boundary: where primes and ideals come from

- Hurwitz--Maass ideals and polarization-only unit transitions come from maximal-order class, narrow-sign, and unit structure. They exist independently of phase 1 and must be tagged `HURWITZ_MAASS_STRUCTURE`.
- Phase 1 representation tests are intended to provide candidate rational primes for unexpected stable torsion subrepresentations. A future implementation must factor each rational prime in `O_F` and retain the actual prime ideal and stable line/plane data. Completeness of the current RM phase-1 proof is still a gap.
- Scalar `L_1,L_2` transitions are node-dependent rational Lagrangians in `B[ell]` or `B[ell^2]`; they cannot be computed once from the initial `A` and reused without proof.
- Nonmaximal conductor primes belong to the excluded ascent/descent interface, not to this slice.

## Classification of the note's named steps

| Named step or assertion | Status | Disposition |
|---|---|---|
| Polarization/isogeny correspondence (`prop:mumford`) | **VERIFIED in scope** | Use Mumford with `beta` and the group scheme over `k` stated explicitly. |
| Reduction to maximal RM (`prop:max-order`) | **OUT_OF_SCOPE / GAP** | “Well known” proof does not establish isotropy, intermediate principal polarizations, or the conductor filtration. |
| Three-stage maximal/nonmaximal reduction corollary | **OUT_OF_SCOPE / GAP** | Depends on the preceding proposition and on unproved descending completeness. |
| Hurwitz--Maass construction from `(a,alpha)` | **VERIFIED in scope** | `A[a]` is base-field rational and maximal isotropic; retain the marked embedding. |
| Hurwitz--Maass finite enumeration | **VERIFIED after correction** | Use `ker(Sq)` extended by positive units, never the malformed later phrase. |
| Cyclic/torsion/inert type definition | **GAP** | “Split = norm p” incorrectly includes ramified primes; ramified local structure needs its own case. |
| Bound `m<=(e-1)/2` in the local normal form | **FALSE AS WRITTEN** | The valid pure-torsion case `e=2,K=A[p]` needs `m=e/2`; use `m<=floor(e/2)` after proof. |
| Inert factorization (`prop:factor-inert`) | **CORRECTABLE / NEEDS_CITATION** | Plausible after a split/inert/ramified partition and explicit pairing proof. |
| Cyclic factorization (`prop:factor-cyclic`) | **GAP** | The displayed subgroup uses the exponent of the cyclic divisor where the torsion exponent is required; it fails already at `e=1`. |
| Cyclic corollary | **GAP** | Depends on the failed proposition; also correct `f o [n]=h o g`. |
| Algorithm step 1: initialize the HM class | **VERIFIED only as arithmetic requests** | Actual target orbit still needs quotient and isomorphism certificates. |
| Algorithm step 2: construct finite `Q` and cyclic neighbors | **GAP** | No finite effective construction or prime-ideal/stable-line bridge is proved; every occurrence must use the current node `B`, not the initial `A`. |
| Algorithm step 3: pass to scalar steps | **SPECIFICATION ONLY** | Safe only after the preceding node closure is effective. |
| Algorithm step 4: compute `L_1,L_2` and recurse | **GAP** | Sets are node-dependent; the RM reducibility proof in `main.tex` still contains `TODO:Proof of finiteness`. |
| Algorithm step 5: output `S=Isog(A)` | **GAP** | No executable equality/deduplication proof and preceding neighbor generators are incomplete. |
| Termination theorem | **GAP** | Faltings finiteness of the answer does not prove finite/effective loops or decidable node equality. |
| Completeness theorem | **GAP** | Blocked by ramified types, cyclic factorization, `Q`, `L_1,L_2`, and nonmaximal reduction. |
| Complex building blocks in `main.tex` | **GAP** | RM preimage, period enumeration, recognition, endomorphism transport, and exact equality are placeholders/heuristics. |

## Paper and implementation disposition

The isolated paper worktree `private/tex/_worktrees/rm-goal/` contains one surviving insertion: the verified fixed-surface polarization proposition above, cited to González--Guàrdia--Rotger. `PAPER_DIFF.patch` is the portable diff. The exact required `latexmk` command passed; `PAPER_LATEXMK.md` preserves the transcript. No Hurwitz--Maass completeness or graph claim was inserted.

The throwaway Sage implementation is `prototype.sage`; its reproducible command, fixtures, exact outputs, independent labels, and validation are in `VALIDATION.md`. The recommended production extraction is specified unequivocally in `NEXT_SLICE.md`.

## Validation disposition

The 13-test public suite and an independent clean process rerun pass.  A
black-box validator checked every positive fundamental discriminant through
`D=260`: 78 fields and 148 reconstructed HM records passed exact ideal-square,
sign, norm, class, and pair-inequivalence checks.  The cardinalities agree
with the elementary genus-theory value `2^(omega(D)-1)` and with
`#Cl+(F)[2]`.  Magma independently reproduced the ordinary and narrow
class-group invariant factors for the seven diagnostic fields.  The exact
canonical ideals and generators were reconstructed through a separate
control path but share Sage/Pari, so `VALIDATION.md` labels that evidence
accordingly rather than calling it an independent CAS result.

The adversarial pass found no arithmetic or paper blocker.  The two formal
review findings were corrected test-first: every HM record now carries its
exact ideal provenance locally, and a nonmaximal rejection has an
order-rejection-only theorem payload with class/unit enumeration marked
`NOT_RUN`.  The independent clean and black-box/Magma agents then reran on the
amended bytes and passed.  Two prototype-only limitations remain: argparse
help is textual rather than result JSON, and the normalization loop has a
defensive rather than proved step cap.  `NEXT_SLICE.md` resolves these at the
production type/API boundary.  No validation in this run constructs or
identifies a quotient PPAS; that layer remains explicitly unvalidated.

## Future-case interfaces

A production caller must distinguish at least:

- `SUPPORTED_EXACT_MAXIMAL_RM`;
- `UNSUPPORTED_NONMAXIMAL_ORDER` with conductor/index data;
- `UNSUPPORTED_RM_NOT_DEFINED_OVER_BASE`;
- `UNSUPPORTED_LARGER_RATIONAL_ENDOMORPHISM_RING`;
- `UNKNOWN_GEOMETRIC_ENDOMORPHISMS` (metadata uncertainty, never silently scoped in);
- `UNKNOWN_ANALYTIC_QUOTIENT`; and
- `UNKNOWN_PPAS_COLLISION`.

This preserves ascent/descent, descent of RM, and larger-endomorphism cases for future work without making an unsupported mathematical claim today.
