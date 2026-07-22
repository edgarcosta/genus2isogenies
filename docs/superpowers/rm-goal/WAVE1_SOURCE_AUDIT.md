# Wave 1 primary-source audit: maximal real-quadratic RM

Date: 2026-07-22

## Verdict

The safe positive slice is narrower than `note.tex`'s advertised full graph
algorithm.  For a geometrically simple PPAS `(A, lambda)/k` together with a
chosen embedding `iota: O_F -> End_k(A)` such that `iota(O_F) = End_k(A)`, the
following are source-supported:

1. the polarization/isogeny correspondence through totally positive symmetric
   endomorphisms and maximal isotropic kernels;
2. enumeration of principal polarizations on the **fixed underlying** `A` by
   `O_F^{x,+}/(O_F^x)^2` (with the action, not merely a count, made explicit);
3. enumeration of Hurwitz--Maass *group elements* by pairs `(p, alpha)`, with
   `p^2 = (alpha)` and `alpha >> 0`, modulo `(p,alpha) ~
   (beta p,beta^2 alpha)`; equivalently, an extension of
   `ker(Sq: Cl(F) -> Cl^+(F))` by `O_F^{x,+}/(O_F^x)^2`.

This does **not** yet prove the full cyclic/torsion/inert BFS in `note.tex`, the
finite construction of `Q,L_1,L_2`, or the paper's analytic reconstruction
interface.  Those remain GAPs below.  In particular, arithmetic
Hurwitz--Maass representatives enumerate a finite group without duplicates,
but their action on a special moduli point can have a stabilizer; target PPASs
must still be deduplicated by polarized `k`-isomorphism.

## Scope and distinctions that must appear in every theorem/API

- `End_k(A)` means endomorphisms defined over the base number field.  The
  scoped equality is an equality **through the chosen embedding**
  `iota: O_F -> End_k(A)`, not an abstract ring-isomorphism.
- `End_{kbar}(A)` may be larger.  No claim below identifies it with `O_F` or
  extends the polarization classification to its extra elements.
- The reference principal polarization `lambda` defines the Rosati involution
  `x^dagger = lambda^{-1} x^vee lambda`.
- In the scoped case, `End_k(A)^dagger = End_k(A) = iota(O_F)`.  Indeed Rosati
  is positive, and the nontrivial automorphism of a real quadratic field is
  not positive (apply it to an element of negative field norm).  Positivity is
  [Mumford, section 21, Theorem 1, pp. 192--193](https://wstein.org/edu/Fall2003/252/references/mumford-abvar/Mumford-Abelian_Varieties.pdf).
- A different chosen RM embedding is different input.  Transport along a
  `k`-isogeny `f` is `x |-> f x f^{-1}` in rational endomorphisms.  The note's
  blanket “unique bijection” is safe only after this embedding/transport is
  stated; it is not a theorem for an unspecified RM subfield inside a larger
  endomorphism algebra.

## Primary sources actually checked

These citations, rather than agreement between local drafts, support the
accepted results.

1. D. Mumford, *Abelian Varieties* (1970):
   - section 20(3), p. 190 identifies `NS(A)_Q` with Rosati-symmetric rational
     endomorphisms;
   - section 21, Theorem 1, pp. 192--193 proves positivity of Rosati;
   - Application III, pp. 208--210 identifies the ample cone with the totally
     positive cone;
   - section 23(1), p. 228 gives functoriality of the commutator pairing;
   - section 23, Corollary to Theorem 2, p. 231 gives descent through an
     isotropic kernel;
   - section 23, Theorem 4, pp. 233--234 gives the order criterion for a
     maximal isotropic subgroup.
   [Primary scan](https://wstein.org/edu/Fall2003/252/references/mumford-abvar/Mumford-Abelian_Varieties.pdf),
   [bibliographic record](https://openlibrary.org/books/OL21166356M).

2. C. Birkenhake and H. Lange, *Complex Abelian Varieties*, 2nd ed. (2004),
   Proposition 5.2.1 and Theorem 5.2.4, pp. 120--121: integral
   Neron--Severi/symmetric-endomorphism correspondence for a principal
   reference polarization and the positive-cone criterion.
   [Chapter DOI](https://doi.org/10.1007/978-3-662-06307-1_7).
   The precise orbit formulation used here is also stated in H. Lange,
   “Principal polarizations on products of elliptic curves,” Proposition 2.1,
   p. 3: principal polarizations correspond to totally positive symmetric
   automorphisms modulo `a |-> u^dagger a u`.
   [arXiv:math/0508458v1](https://arxiv.org/pdf/math/0508458v1).

3. G. van der Geer, *Hilbert Modular Surfaces* (1988):
   - section I.4, pp. 12--13 defines
     `Sq: Cl(F) -> Cl^+(F), [p] |-> [p^2]`, chooses ideal representatives of
     its kernel, and gives the Hurwitz--Maass quotient and its unit kernel;
   - section IX.1, Definition 1.5 and Proposition 1.6, pp. 208--209 distinguishes
     an ordered polarization module from a fixed polarization and notes that
     positive scaling changes the polarization.
   [Book DOI](https://doi.org/10.1007/978-3-642-61553-5),
   [section I.4 chapter](https://link.springer.com/book/10.1007/978-3-642-61553-5).

4. G. Faltings, “Endlichkeitssatze fur abelsche Varietaten uber
   Zahlkorper,” *Invent. Math.* 73 (1983), 349--366,
   [DOI](https://doi.org/10.1007/BF01388432).  This is a valid source for
   finiteness consequences, but it does not repair the missing effective
   finiteness/computability arguments for `Q,L_1,L_2` in the note.

## Local note versus arXiv v1 (one comparison; non-independent copies)

`private/tex/note.tex` was compared byte-for-byte, after decompressing the
official TeX source, with
[arXiv:2405.19820v1](https://arxiv.org/abs/2405.19820v1) (submitted 30 May
2024).  There are no mathematical differences.  The only differences are:

- the local abstract calls the types “elementary” in quotation marks and says
  the note is intended for a future joint paper;
- arXiv v1 instead says it is ongoing joint work and not intended for separate
  publication;
- the local bibliography adds that future joint paper.

The two texts are versions of the same work and are **not independent
validation** of any claim.

## Dependency and claim audit

Status meanings: VERIFIED = checked against the cited primary result and the
scoped hypotheses; NEEDS_CITATION = plausible but the manuscript lacks an
adequate exact citation/proof; GAP = false as written, proof incomplete, or an
algorithmic implication is missing; OUT_OF_SCOPE = explicitly excluded from
this implementation slice.

| Claim/dependency | Status | Audit result |
|---|---|---|
| Isogenies to PPAVs correspond to `(beta,G)` with `beta` symmetric totally positive and `G` maximal isotropic in `A[beta]` | VERIFIED | Mumford section 23, pp. 228, 231, 233--234, plus the positive cone on pp. 208--210.  Rationality requires `beta` and `G` over `k`; the note's Proposition 1 must say this explicitly. |
| `deg(beta)=deg(f)^2` | VERIFIED | Take degrees in `f^vee lambda' f = lambda beta`; both principal polarizations have degree one. |
| Rosati fixes the scoped `O_F` | VERIFIED | Uses `End_k(A)=iota(O_F)` and Mumford section 21, Theorem 1.  This does not assert anything about extra geometric endomorphisms. |
| Every principal polarization on fixed `A/k`, up to polarized `k`-isomorphism, is represented by a totally positive unit modulo the automorphism action | VERIFIED | Lange Proposition 2.1.  In this scope `Aut_k(A)=O_F^x`, dagger is identity, so the set is `O_F^{x,+}/(O_F^x)^2`. |
| Changing the reference polarization | VERIFIED | If `lambda_e = lambda o e`, then `x^{dagger_e}=e^{-1}x^dagger e`.  On the commutative `End_k(A)=O_F` this remains identity; on a larger `End_{kbar}(A)` it need not. |
| `A[p]` is `k`-rational for an invertible `O_F`-ideal `p` | VERIFIED | It is the intersection of kernels of base-field endomorphisms in `p`. |
| If `p^2=(alpha)`, `alpha >> 0`, then `A[p]` is maximal isotropic in `A[alpha]` | VERIFIED | After `k -> C`, the `O_F`-stable homology lattice is projective of rank two.  Self-adjointness gives integrality of the Riemann form on `p Lambda x p^{-1} Lambda`; `#A[p]=N(p)^2=sqrt(#A[alpha])`.  Apply Mumford section 23, Theorem 4.  This short lattice proof should be included or cited in the paper; van der Geer I.4 alone is too terse for the algebraic assertion. |
| The Hurwitz--Maass quotient still has `End_k=O_F` | VERIFIED | The `O_F` action descends through the stable kernel.  Rational endomorphism algebras are transported by the `k`-isogeny, and the resulting integral order contains the already maximal `O_F`. |
| Correct finite Hurwitz--Maass datum | VERIFIED | van der Geer I.4, pp. 12--13 gives an exact extension with quotient `ker(Sq:Cl(F)->Cl^+(F))` and kernel `O_F^{x,+}/(O_F^x)^2`; equivalently use pair classes `(p,alpha)/(beta p,beta^2 alpha)`. |
| Note's later phrase “classes in `Cl^+(F)` that are squares of elements of `Cl(F)`” | GAP | It is ill-typed/wrong.  Replace it by “classes `[p] in ker(Sq:Cl(F)->Cl^+(F))`,” then multiply a fixed positive generator by representatives of `O_F^{x,+}/(O_F^x)^2`. |
| Initial condition `[p]^2=1 in Cl^+(F)` versus later enumeration | VERIFIED after correction | For an individual ideal it is exactly `p^2=(alpha)` with `alpha>>0`.  Quotienting pairs by arbitrary principal rescaling produces the extension above.  It is not legitimate simply to identify the data with the malformed later phrase (nor canonically with `Cl^+(F)[2]`, although van der Geer proves the cardinalities agree). |
| Arithmetic HM representatives give the orbit without duplicate PPASs | GAP | They give distinct elements of `Gamma_m/Gamma`; a point can have a nontrivial stabilizer.  The analytic engine must return a polarized-isomorphism key and the caller must deduplicate. |
| Reduction from a nonmaximal order to maximal RM (note Proposition 2.1) | OUT_OF_SCOPE | The proof is labelled “well known” and omits the lattice/polarization justification.  It must not be used by the maximal-only prototype. |
| RM not defined over `k`, or an unspecified RM subfield of a larger `End_k^0` | OUT_OF_SCOPE | Kernel rationality and uniqueness of transported RM do not follow from the scoped argument. |
| QM, `M_2(Q)`, products, and larger rational endomorphism rings | OUT_OF_SCOPE | The fixed-surface polarization quotient is then a noncommutative congruence problem, not a unit-square quotient. |
| Local cyclic/torsion normal form in the full classification | GAP | The stated bound `m(p) <= (e(p)-1)/2` excludes the valid pure-torsion case `e=2, K=A[p]`, which requires `m=e/2`. |
| Proposition `factor-cyclic` | GAP | Its proof defines `e''` as the exponent in `q`; the subgroup `A[p^{e''}]` requires the torsion divisor `r`.  The older draft had `r`, confirming a material typo.  The finite-group calculation is otherwise only asserted. |
| Factorization of the inert type | NEEDS_CITATION | Plausible from the primary decomposition, but the induced pairing and maximality are not sourced in the note. |
| Existence/effective computation of finite `Q` in Algorithm 5.4 step 2 | GAP | No proof that the required divisibility basis is finite or computable is supplied.  “Will rely on Dieulefait” is not an algorithm or termination proof. |
| Effective finite sets `L_1,L_2` | GAP | `private/tex/main.tex` itself leaves the RM reducibility proof's finiteness as `TODO`.  This cannot be accepted into a paper-ready theorem. |
| Termination/correctness of full Algorithm 5.4 | GAP | Faltings finiteness of the state space does not prove that the unspecified neighbor generators are finite/effective; it also does not repair the two local factorization defects above. |
| Complex building blocks / quotient recognition in `private/tex/main.tex` | GAP | The RM preimage, period enumeration, rational recognition, endomorphism transport, and certified equality tests are placeholders/questions, not proved algorithms. |
| Phase-1 prime list versus RM structural ideals | NEEDS_CITATION | Phase 1 should supply only exceptional reducible Galois primes/ideals.  Hurwitz--Maass ideal classes and positive-unit classes are unconditional RM structure and must not be conflated with that list.  The exact phase-1 completeness theorem is still missing. |

## Exact fixed-surface polarization statement ready for the paper

Let `(A,lambda)/k` satisfy the scoped hypotheses.  For every `e in
O_F^{x,+}`, `lambda_e = lambda o iota(e)` is a principal polarization over
`k`.  Every principal polarization over `k` is uniquely of this form before
quotienting by automorphisms.  Two such polarizations `lambda_e` and
`lambda_e'` are `k`-isomorphic exactly when

`e' = u^dagger e u = u^2 e` for some `u in O_F^x`.

Hence

`PP_k(A) / ~=  O_F^{x,+} / (O_F^x)^2`.

For a real quadratic field this quotient is trivial when a unit of negative
norm exists, and otherwise has order two.  This is a statement about the fixed
underlying `A` and `k`-defined polarizations; it is not a statement about
geometric polarizations arising from `End_{kbar}(A)`.

## Minimum contract for the safe prototype/analytic boundary

The exact-arithmetic prototype may safely accept only `(F,O_F)` and reject a
nonmaximal order.  It can return canonical representatives of:

- `ker(Sq: Cl(F) -> Cl^+(F))`;
- `O_F^{x,+}/(O_F^x)^2`;
- pair classes `(p,alpha)` with machine-checked `p^2=(alpha)` and both real
  embeddings of `alpha` positive.

To continue a geometric BFS, an analytic quotient engine must additionally
return (and certify): the target PPAS and principal polarization, the
transported chosen `O_F` embedding, a polarized `k`-isomorphism/deduplication
key, and an explicit `unknown/unsupported` result when rational recognition or
endomorphism certification fails.  None of these analytic certificates is
provided by the cited Hurwitz--Maass group theorem itself.
