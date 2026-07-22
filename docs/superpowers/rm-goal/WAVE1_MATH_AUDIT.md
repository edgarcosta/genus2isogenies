# Wave 1 mathematical audit: maximal real-quadratic RM

## Verdict

The safe frozen slice is **not** the full algorithm in `note.tex`. It is the
finite arithmetic of Hurwitz--Maass pairs, together with the classification of
principal polarizations on a fixed underlying surface. Both admit short proofs
under the exact scoped hypothesis

\[
 (A,\lambda,\iota),\qquad
 \iota:\mathcal O_F\xrightarrow{\sim}\operatorname{End}_k(A),
\]

where (A/k) is geometrically simple, (F) is real quadratic, and
(\mathcal O_F) is maximal. The embedding (\iota), the polarization
(\lambda), and the base field are part of the data.

The claimed full enumeration theorem (`private/tex/note.tex:471-490`) is **not
paper-ready**. It has independent blockers: ramified primes are misclassified,
one local normal-form bound and one kernel formula are wrong, and the finiteness
and computability of (Q,L_1,L_2) are not proved.

Status words below mean:

- **VERIFIED**: a valid argument is given here under the frozen hypotheses.
- **CORRECTABLE**: the intended claim appears valid after the stated repair.
- **GAP**: the current proof cannot support the claim.
- **OUT_OF_SCOPE**: no positive claim should be made in this run.

## 1. Scope, endomorphisms, and Rosati

### RM-01 — the relevant ring is (\operatorname{End}_k(A)), not the geometric ring — VERIFIED

The note consistently begins with objects defined over (k)
(`note.tex:72-84`), but later material alternates between
(\operatorname{End}_k(A)) and (\operatorname{End}_{\bar k}(A)). The frozen
input must assert an exact chosen isomorphism
(\iota:\mathcal O_F\simeq\operatorname{End}_k(A)). It must **not** infer
(\operatorname{End}_{\bar k}(A)=\mathcal O_F).

Geometric extra endomorphisms may exist without being usable over (k). The
proofs below use only (k)-endomorphisms and (k)-isomorphisms. Any algorithm
that deduplicates by complex/absolute Igusa invariants alone would therefore
collapse descent forms and is unsafe; the typical implementation's separate
twist-pinning step is an interface oracle, not a proof that absolute invariants
are a (k)-isomorphism key.

Cases with (\operatorname{End}_k(A)\supsetneq\iota(\mathcal O_F)), RM not
defined over (k), nonmaximal orders, products, (M_2(\mathbb Q)), QM, and CM
classification are **OUT_OF_SCOPE**. Geometric extra endomorphisms are opaque
metadata: they confer no extra allowed operations in this slice.

### RM-02 — the chosen (\mathcal O_F) is Rosati-fixed — VERIFIED

The Rosati involution associated with any principal polarization restricts to a
positive involution of
(\operatorname{End}_k^0(A)=F). It cannot be the nontrivial automorphism of a
real quadratic field: for an (x\in F) with (N_{F/\mathbb Q}(x)<0), positivity
would give
(\operatorname{Tr}_{F/\mathbb Q}(x x^\dagger)=2N(x)<0). Hence it is the
identity. Thus

\[
 \operatorname{End}_k(A)^\dagger=\operatorname{End}_k(A)=\iota(\mathcal O_F).
\]

This justifies replacing the note's weaker and potentially misleading
(F=\operatorname{End}_k(A)^\dagger\otimes\mathbb Q)
(`note.tex:156-169`) by the exact scoped hypothesis. The suggested extension to
QM and matrix algebras (`note.tex:171-177`) is **OUT_OF_SCOPE**: in a
noncommutative algebra the Rosati-fixed set is generally not a subring.

### RM-03 — transported embeddings must be retained — VERIFIED

For an isogeny \(f:A\to B\) between scoped maximal-RM surfaces, transport gives
\(x\mapsto fxf^{-1}\) on rational endomorphism algebras. Since both integral
endomorphism rings are the unique maximal order of (F), this transports
(\mathcal O_F) to (\operatorname{End}_k(B)). Consequently (f) is
(\mathcal O_F)-linear for the **transported** embedding and its kernel is
(\mathcal O_F)-stable. This is a conclusion, not an unrecorded assumption as
at `note.tex:274-283`.

Conjugating the initial abstract embedding by the nontrivial automorphism of
(F) conjugates every ideal label. An implementation must either retain the
chosen embedding or explicitly quotient its output by this involution.

## 2. Principal polarizations on the fixed surface

### POL-01 — exact classification — VERIFIED

Fix the scoped ((A,\lambda,\iota)), and let
(U=\mathcal O_F^\times), (U^+=\mathcal O_F^{\times,+}). Then the principal
polarizations on the fixed underlying (A), defined over (k), up to
polarized (k)-isomorphism are in bijection with

\[
                         U^+/U^2.                 \tag{1}
\]

Proof. Relative to (\lambda), a second polarization is
(\lambda_\epsilon=\lambda\circ\iota(\epsilon)) for a Rosati-symmetric
totally positive endomorphism (\epsilon). It is principal exactly when
(\deg(\epsilon)=N_{F/\mathbb Q}(\epsilon)^2=1), hence exactly when
(\epsilon\in U^+). Every (k)-automorphism is a unit (u\in U), and

\[
 u^*\lambda_\epsilon
 =u^\vee\lambda\epsilon u
 =\lambda\,u^\dagger\epsilon u
 =\lambda\,\epsilon u^2.
\]

Therefore the orbit relation is multiplication by (U^2), proving (1). This
is the scoped specialization of the standard polarization/positive-symmetric
endomorphism correspondence (Mumford; Birkenhake--Lange, *Complex Abelian
Varieties*, Remark 5.2.5).

For a real quadratic field this quotient has order one when a unit of negative
norm exists, and order two otherwise.

### POL-02 — changing the reference polarization — VERIFIED

For (\lambda_\epsilon=\lambda\epsilon),

\[
 x^{\dagger_\epsilon}=\epsilon^{-1}x^\dagger\epsilon.
\]

On (\operatorname{End}_k^0(A)=F) this is still the identity because (F) is
commutative. It may change the involution and fixed subspace inside a larger
(\operatorname{End}_{\bar k}^0(A)); classifying that geometric effect is
**OUT_OF_SCOPE**. The paper must not say merely that “Rosati depends on
(\lambda)” without recording this scoped invariance.

## 3. Hurwitz--Maass isogenies

### HM-01 — construction from a pair — VERIFIED

Define a Hurwitz--Maass pair to be

\[
 (\mathfrak a,\alpha),\qquad
 \mathfrak a^2=(\alpha),\quad \alpha\gg0,
\]

with (\mathfrak a) an invertible fractional (\mathcal O_F)-ideal. Use an
integral representative to define
(A[\mathfrak a]=\bigcap_{x\in\mathfrak a}\ker\iota(x)).
Because the action is over (k), this group scheme is (k)-rational. Locally,
(T_\ell A) is free of rank two over
(\mathcal O_F\otimes\mathbb Z_\ell); the relation
(\mathfrak a^2=(\alpha)), self-adjointness of (\mathcal O_F), and the
Riemann form show that (A[\mathfrak a]) is isotropic in (A[\alpha]). Its
order is (N(\mathfrak a)^2=\sqrt{\#A[\alpha]}), so it is maximal isotropic.
The standard Mumford quotient lemma (`note.tex:101-124`) therefore gives a
unique principal polarization on (A/A[\mathfrak a]).

The target still has exact (k)-endomorphism ring (\mathcal O_F): transported
(\mathcal O_F) acts on it, and a (k)-isogenous variety has the same rational
endomorphism algebra (F), so maximality leaves no larger order in (F).

This is also the ideal condition in the classical Hurwitz--Maaß extension; see
Hirzebruch--Zagier, §III.3.1, [*Classification of Hilbert Modular
Surfaces*](https://hirzebruch.mpim-bonn.mpg.de/100/1/61_Classification%20of%20Hilbert%20modular%20spaces.pdf).
The typo at `note.tex:235` must read (\lambda_B), not (\lambda_A).

### GAP-HM-ENUM-01 — kernel of squaring is not image of squaring

The initial condition at `note.tex:228-233` is correct. The enumeration at
`note.tex:248-255` is not equivalent to it.

There is a well-defined homomorphism

\[
 \operatorname{Sq}:\operatorname{Cl}(F)\longrightarrow
 \operatorname{Cl}^+(F),\qquad
 [\mathfrak a]\longmapsto[\mathfrak a^2]^+ .       \tag{2}
\]

It is well-defined precisely because changing \(\mathfrak a\) by
\((x)\mathfrak a\) changes its square by the *narrow-principal* ideal
\((x^2)\). Eligible ordinary ideal classes are
(\ker(\operatorname{Sq})). They are not narrow classes in
(\operatorname{im}(\operatorname{Sq})), nor are they simply
(\operatorname{Cl}^+(F)[2]).

Two concrete counterexamples:

- (F=\mathbb Q(\sqrt{10})):
  (\operatorname{Cl}(F)\simeq\operatorname{Cl}^+(F)\simeq C_2) and a
  negative-norm unit exists. Thus (2) has kernel (C_2) and trivial image.
  Enumerating the image misses the nontrivial eligible ideal class.
- (F=\mathbb Q(\sqrt{229})):
  (\operatorname{Cl}(F)\simeq\operatorname{Cl}^+(F)\simeq C_3), again with a
  negative-norm unit. Squaring is an automorphism, so its kernel is trivial and
  its image is all of (C_3). Enumerating the image adds two ineligible
  classes.

These class-group and unit assertions were also checked directly in Sage.

### HM-02 — the correct finite datum — VERIFIED

Let

\[
 \mathcal H_F=
 \{(\mathfrak a,\alpha):\mathfrak a^2=(\alpha),\ \alpha\gg0\}/\sim,
\]

where
((\mathfrak a,\alpha)\sim(x\mathfrak a,x^2\alpha)) for (x\in F^\times).
Then there is an exact sequence of finite elementary (2)-groups

\[
 1\longrightarrow U^+/U^2\longrightarrow\mathcal H_F
 \longrightarrow\ker(\operatorname{Sq})\longrightarrow1.       \tag{3}
\]

Proof. Projection to the ordinary class of (\mathfrak a) has the displayed
image. Its kernel can be scaled to ((\mathcal O_F,\epsilon)), with
(\epsilon\in U^+), and two such pairs differ by a unit square. Squaring any
pair and scaling by (\alpha^{-1}) gives the identity, so every element has
order at most two.

For duplicate-free enumeration, choose one deterministic ordinary ideal
representative (\mathfrak a_c) for every
(c\in\ker(\operatorname{Sq})), one verified totally positive generator
(\alpha_c) of (\mathfrak a_c^2), and emit

\[
       (\mathfrak a_c,\epsilon\alpha_c),
       \qquad \epsilon\in U^+/U^2.                  \tag{4}
\]

No splitting of (3) is claimed; (4) is a set of representatives. Under
(\operatorname{End}_k(A)=\mathcal O_F), two data in (4) give polarized
(k)-isomorphic quotients only if they are equivalent above: a polarized
isomorphism between the quotients induces multiplication by an (x\in F^\times)
on their RM lattices, forcing
(\mathfrak b=x\mathfrak a) and (\beta=x^2\alpha). Geometric fixed points
coming from extra endomorphisms are irrelevant to this (k)-isomorphism claim.

The subgroup (U^+/U^2) in (3) is exactly the fixed-surface polarization list
from (1); a pair ((\mathcal O_F,\epsilon)) has zero kernel and changes
(\lambda) to (\lambda\epsilon).

## 4. Cyclic, torsion, and inert types

### GAP-TYPE-RAMIFIED-01 — “split” wrongly includes ramified primes

At `note.tex:287-293`, “split” is defined by (N(\mathfrak p)=p). This includes
ramified prime ideals. The asserted model

\[
 A[\mathfrak p^e]\simeq(\mathbb Z/p^e\mathbb Z)^2
\]

is valid for an unramified split prime but false when (p) ramifies. For
example, if (\mathfrak p^2=(p)), then
(A[\mathfrak p^2]=A[p]\simeq(\mathbb Z/p\mathbb Z)^4), not
((\mathbb Z/p^2\mathbb Z)^2). Ramified primes are also absent from the stated
split/inert dichotomy. They require a separate DVR/uniformizer argument before
the type classification can be accepted.

### GAP-TYPE-NORMALFORM-02 — the bound on (m(\mathfrak p)) is false

Even at a genuinely split prime, `note.tex:302` must allow

\[
 0\le m(\mathfrak p)\le\lfloor e(\mathfrak p)/2\rfloor,
\]

not (m\le(e-1)/2). For (e=2), the subgroup
(K=A[\mathfrak p]\subset A[\mathfrak p^2]) is maximal isotropic and has
((n_0,n_1)=(2,0)), hence (m=1=e/2). The intended unramified local normal
form is (K\simeq\mathbb Z/p^{e-m}\oplus\mathbb Z/p^m).

### GAP-FACTOR-CYCLIC-03 — the displayed kernel is wrong

In `note.tex:378-390`, write
(a=v_{\mathfrak p}(\mathfrak q)) and
(b=v_{\mathfrak p}(\mathfrak r)). The term added should be
(A[\mathfrak p^b]), not (A[\mathfrak p^a]). For a genuinely split prime the
candidate repair is

\[
 K'_{\mathfrak p}=A[\mathfrak p^b]
       +p^{,e-m-a-b}K_{\mathfrak p}.               \tag{5}
\]

The current formula already fails for (e=1,m=0,
\mathfrak q=\mathfrak p,\mathfrak r=1): it returns all of
(A[\mathfrak p]), of order (p^2), instead of the cyclic Lagrangian line
(K_{\mathfrak p}), of order (p). Thus the proposition may be repairable,
but its proof is presently invalid. Formula (5) still needs a written local
proof and a ramified analogue.

### RM-TYPE-04 — inert separation — CORRECTABLE

For an inert rational prime (p), the prime ideal is ((p)), so the inert ideal
has a positive rational generator. Orthogonal CRT decomposition then separates
the inert kernel as intended in `note.tex:331-355`. This claim is usable after
the prime categories are changed to **split / inert / ramified**, transported
RM is made explicit, and ideal factorizations are written
((\alpha)=\prod\mathfrak p^{e(\mathfrak p)}).

The cyclic corollary also has a composition-order typo: the diagram gives
(f\circ[n]=h\circ g), not (g\circ h) (`note.tex:399-417` and `:484`). It
depends on GAP-FACTOR-CYCLIC-03 and is not yet paper-ready.

## 5. (Q,L_1,L_2), termination, and completeness

### GAP-Q-01 — (Q) is a specification, not a computed finite set

At `note.tex:429-447`, all occurrences of (A[\mathfrak q]) and quotients of
(A) must refer to the **current node (B)**. “Maximal cyclic subgroup” must be
defined as an (\mathcal O_F/\mathfrak q)-cyclic Lagrangian, including its
prime-power behavior. The note neither proves that a finite divisibility basis
(Q_B) is computable nor connects rational-prime output to individual prime
ideals and stable lines.

A possible future repair is: obtain a finite set of candidate rational primes
from phase 1; factor them in (\mathcal O_F); test the residue-degree-one ideal
representations for (k)-stable lines; then enumerate the finitely many minimal
products whose narrow class is a square. This is not proved in the current
manuscript.

### GAP-L12-02 — (L_1,L_2) are node-dependent and their finiteness is unproved

At `note.tex:450-463`, (A[\ell]) must again be (B[\ell]). The definitions
describe actual scalar 1-step and 2-step kernels, but do not give finite
candidate lists. The RM reducibility proposition in `main.tex:606-656` ends
with the explicit marker `TODO:Proof of finiteness` (`main.tex:644`). Its input
also conflates the RM field with the endomorphism/base field. Consequently it
cannot yet certify (Q_B,L_{1,B},L_{2,B}).

### GAP-TERMINATION-03 — finiteness of the answer does not prove executable termination

Faltings finiteness plus the finite quotient (1) shows that the scoped set of
polarized (k)-isomorphism classes is finite. It does not prove that loops over
undefined (Q,L_1,L_2) terminate, that equality is decidable, or that every
node has finitely many attempted edges. The argument at `note.tex:475-477` is
therefore insufficient.

### GAP-COMPLETENESS-04 — the final theorem is blocked

The induction on cyclic type is a plausible architecture: after inert
separation, remove a nontrivial divisor of the cyclic ideal and finish in the
Hurwitz--Maass orbit. But it depends on GAP-TYPE-RAMIFIED-01,
GAP-FACTOR-CYCLIC-03, GAP-Q-01, and GAP-L12-02. Reduction from/to nonmaximal RM
is additionally **OUT_OF_SCOPE**. No full completeness proposition should be
inserted into the paper in this wave.

The reduction-to-maximal-order proposition (`note.tex:183-218`) is itself a
**GAP/OUT_OF_SCOPE** item: the short “well known” proof does not establish the
claimed isotropy, induced maximal action, intermediate principal polarizations,
or conductor-prime filtration. It must not be used by the frozen slice.

## 6. Phase boundaries and analytic return contract

The candidate sources must remain distinct:

- **Hurwitz--Maass structure:** (\ker(\operatorname{Sq})) and (U^+/U^2).
  These kernels are rational because the full (\mathcal O_F)-action is over
  (k). They must not be filtered through `reducible_ell`.
- **Phase 1 representation tests:** candidate rational primes for noncanonical
  stable lines/planes; after factoring (\ell\mathcal O_F), the algorithm must
  identify the relevant prime ideal and subrepresentation. Bad and ramified
  primes require explicit handling. This phase is presently a GAP.
- **Scalar (L_1,L_2) steps:** actual (k)-rational Lagrangians in
  (B[\ell]) or (B[\ell^2]), recomputed or safely transported per node.

For RM breadth-first enumeration to continue, an analytic quotient engine must
return more than modular invariants. Its exact success value must contain

1. a (k)-model/representation of the target PPAS (B') and its principal
   polarization (\lambda_{B'});
2. the induced chosen embedding
   (\iota_{B'}:\mathcal O_F\simeq\operatorname{End}_k(B')), with exact action
   matrices on differentials/homology sufficient for the next ideal kernel;
3. an exact kernel/degree/ideal certificate and the pullback identity
   (f^*\lambda_{B'}=\lambda_B\alpha);
4. a (k)-polarized-isomorphism certificate/key for deduplication, including
   twist/descent data rather than only complex Igusa invariants; and
5. an explicit `unsupported`/`unknown` result when exact recognition or the RM
   action cannot be certified.

The complex-building-block text (`main.tex:716-742`) is currently a placeholder
with heuristic rounding and an unchecked representation factor, so it does not
satisfy this contract.

## 7. Frozen recommendation

Implement only the following exact arithmetic slice in Wave 2:

1. accept a real quadratic maximal order (\mathcal O_F); explicitly reject a
   nonmaximal order;
2. compute (U^+/U^2);
3. compute (2) and enumerate (\ker(\operatorname{Sq})), **not its image**;
4. emit deterministic representatives (4), each with machine-checked
   (\mathfrak a^2=(\alpha)) and positivity under both real embeddings; and
5. expose the polarization representatives (1) separately while recording
   that they are the unit-kernel subgroup of (3).

The prototype may describe the quotient request
((A,\lambda,\iota;\mathfrak a,\alpha)), but must not claim to compute the full
RM graph or the quotient surface unless the analytic return contract above is
met. Recommended adversarial fixtures are
(\mathbb Q(\sqrt3)) (unit contribution but trivial ordinary class group),
(\mathbb Q(\sqrt{10})) (nontrivial kernel of Sq but trivial image), and
(\mathbb Q(\sqrt{229})) (trivial kernel but surjective Sq).

Paper-ready scoped propositions are: RM-02, POL-01/POL-02, HM-01, and HM-02,
after primary citations and notation cleanup. Every broader algorithmic claim
listed as a GAP above must remain out of the insertion.
