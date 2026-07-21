// IsogenyPrimes.m
//
// Isogeny primes / congruence primes engine for elliptic curves over Q and
// number fields. See docs/specs/billerey-engine-design.md
// for the denotational semantics this file implements.

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

//////////////////////////////////////////////////////////////////////////
// Internal helpers (plain functions, not intrinsics)
//////////////////////////////////////////////////////////////////////////

function CharpolyFromTraceAndNorm(a, n)
    Zx<t> := PolynomialRing(Integers());
    return t^2 - a*t + n;
end function;

// A global integral model of E, minimal where Magma can certify one (over Q
// always; over a number field only when the class number allows it).
// MinimalModel(E) requires class number 1 over a genuine number field and
// raises a require-failure otherwise, so we fall back to IntegralModel.
function NormalizedModel(E)
    try
        Enorm := MinimalModel(E);
        return Enorm;
    catch err
        return IntegralModel(E);
    end try;
end function;

// Prime support of the conductor of E, which MUST already be a normalized
// (global integral) model: rational primes over Q, prime ideals over a
// number field. Never call this on the curve as submitted by the caller.
function ConductorSupport(E)
    N := Conductor(E);
    if Type(BaseRing(E)) eq FldRat then
        return PrimeDivisors(N);
    end if;
    return [ pr[1] : pr in Factorization(N) ];
end function;

//////////////////////////////////////////////////////////////////////////
// FrobeniusCharpoly
//
// Magma does not allow user intrinsics to declare a parameter of type Any
// (verified: "Illegal intrinsic"), so the P::Any of the design spec is
// realized as two overloads, one per admissible type of P.
//////////////////////////////////////////////////////////////////////////

intrinsic FrobeniusCharpoly(E::CrvEll, P::RngIntElt) -> RngUPolElt
{The characteristic polynomial x^2 - a*x + P of Frobenius at the good
rational prime P, for E defined over the rationals.}
    require Type(BaseRing(E)) eq FldRat :
        "FrobeniusCharpoly: E must be defined over the rationals when P is a rational prime";
    require IsPrime(P) : "FrobeniusCharpoly: P must be a prime integer";
    require Valuation(Discriminant(E), P) eq 0 :
        "FrobeniusCharpoly: E must have good reduction at P";
    Ered := ChangeRing(E, GF(P));
    a := TraceOfFrobenius(Ered);
    return CharpolyFromTraceAndNorm(a, P);
end intrinsic;

intrinsic FrobeniusCharpoly(E::CrvEll, P::RngOrdIdl) -> RngUPolElt
{The characteristic polynomial x^2 - a*x + Norm(P) of Frobenius at the good
prime ideal P, for E defined over a number field.}
    require Type(BaseRing(E)) eq FldNum :
        "FrobeniusCharpoly: E must be defined over a number field when P is a prime ideal";
    require IsPrime(P) : "FrobeniusCharpoly: P must be a prime ideal";
    // E need not be locally minimal at P (class number > 1 fields have no
    // global minimal model); MinimalModel(E,P) always succeeds for a single
    // prime and makes the discriminant valuation the true good-reduction test.
    Emin := MinimalModel(E, P);
    require Valuation(Discriminant(Emin), P) eq 0 :
        "FrobeniusCharpoly: E must have good reduction at P";
    Ered := Reduction(Emin, P);
    a := TraceOfFrobenius(Ered);
    return CharpolyFromTraceAndNorm(a, Norm(P));
end intrinsic;

//////////////////////////////////////////////////////////////////////////
// HasPrimeIsogeny: the pinned shared validated-kernel isogeny primitive.
//
// IsogeniesPrimeDegree and Isogenies do not exist on Magma 2.29-8, and
// IsogenousCurves/IsIsogenous error on CrvEll over FldNum; IsogenyFromKernel
// is the one primitive verified to work over both Q and number fields, so
// it is the validator for hand-built kernel-polynomial candidates below.
//////////////////////////////////////////////////////////////////////////

// Index-subsets of degs (a SeqEnum of positive integers, in bijection with
// a division polynomial's K-irreducible factors -- each already a
// Galois-stable root set, so any union of factors is Galois-stable too)
// summing exactly to target. Suffix-sum pruning means an unreachable target
// costs O(#degs) and enumerates nothing, which is what keeps the large
// Mazur primes (67, 163) cheap on curves with no isogeny of that degree.
function IndexSubsetsSummingTo(degs, target)
    n := #degs;
    suffix := [0 : i in [1 .. n+1]];
    for i in [n .. 1 by -1] do
        suffix[i] := suffix[i+1] + degs[i];
    end for;
    if target eq 0 then
        return [ [Integers()|] ];
    end if;
    if suffix[1] lt target then
        return [ Parent([Integers()|]) | ];
    end if;

    function Rec(i, remaining)
        if remaining eq 0 then
            return [ [Integers()|] ];
        end if;
        if i gt n or suffix[i] lt remaining then
            return [ Parent([Integers()|]) | ];
        end if;
        results := Rec(i + 1, remaining);            // skip factor i
        if degs[i] le remaining then
            withI := Rec(i + 1, remaining - degs[i]);
            results cat:= [ [i] cat s : s in withI ];  // take factor i
        end if;
        return results;
    end function;

    return Rec(1, target);
end function;

// A handful of good-reduction Frobenius charpolies of E (integer
// coefficients, via FrobeniusCharpoly), skipping the rational prime skip.
// Works uniformly over FldRat and FldNum (including degree-one FldNum, e.g.
// RationalsAsNumberField()): NormalizedModel is not assumed here, so callers
// pass in whatever model they already have.
function SampleFrobeniusCharpolies(E, count, skip)
    K := BaseRing(E);
    Disc := Discriminant(E);
    data := [];
    ell := 2;
    tried := 0;
    while #data lt count and tried lt 2000 do
        if IsPrime(ell) and ell ne skip then
            tried +:= 1;
            if Type(K) eq FldRat then
                if Valuation(Disc, ell) eq 0 then
                    Append(~data, FrobeniusCharpoly(E, ell));
                end if;
            else
                for pr in Decomposition(Integers(K), ell) do
                    P := pr[1];
                    if Valuation(Disc, P) eq 0 then
                        Append(~data, FrobeniusCharpoly(E, P));
                    end if;
                end for;
            end if;
        end if;
        ell := NextPrime(ell);
    end while;
    return data;
end function;

// Necessary-condition short circuit: if E[p] were reducible, every good
// Frobenius would have a mod-p eigenvalue on the stable line, so its
// charpoly would be reducible mod p. One good charpoly that is IRREDUCIBLE
// mod p already proves p is not an isogeny prime, without ever forming the
// p-division polynomial. Same style of test as this repo's own
// genus2isogenies.py reducible_ell and the design spec's branch-2 Frobenius
// filter; purely a performance short-circuit (the 67- and 163-division
// polynomials are otherwise minutes-scale to factor on some curves -- see
// the design spec's Risks section on capping pathological enumeration).
function FrobeniusExcludesPrime(E, p)
    Fp := GF(p);
    Fpx := PolynomialRing(Fp);
    for cp in SampleFrobeniusCharpolies(E, 12, p) do
        if IsIrreducible(Fpx ! cp) then
            return true;
        end if;
    end for;
    return false;
end function;

intrinsic HasPrimeIsogeny(E::CrvEll, p::RngIntElt) -> BoolElt, SeqEnum
{Whether E has an isogeny of prime degree p defined over its base field, and
the full sequence of validated codomains (one per validated candidate
kernel polynomial; may have more than one entry, e.g. two independent
p-isogenies out of E). Candidates: for p = 2, the linear factors of the
2-division polynomial (rational 2-torsion x-coordinates); for odd p,
products of Galois-stable factor subsets of the p-division polynomial of
total degree (p-1)/2. A candidate is validated iff
IsogenyFromKernel(E, candidate) succeeds.}
    require IsPrime(p) : "HasPrimeIsogeny: p must be prime";

    if FrobeniusExcludesPrime(E, p) then
        return false, [];
    end if;

    f := DivisionPolynomial(E, p);
    factors := [ g[1] : g in Factorization(f) ];

    if p eq 2 then
        candidates := [ h : h in factors | Degree(h) eq 1 ];
    else
        target := (p - 1) div 2;
        subsets := IndexSubsetsSummingTo([ Degree(h) : h in factors ], target);
        error if #subsets gt 10000, Sprintf(
            "HasPrimeIsogeny: candidate kernel enumeration exceeded its safety cap at p = %o " *
            "(%o Galois-stable subsets); this looks like the pathological case flagged in the " *
            "design spec's Risks section.", p, #subsets);
        candidates := [ &*[ factors[i] : i in s ] : s in subsets ];
    end if;

    codomains := [];
    for cand in candidates do
        try
            Append(~codomains, IsogenyFromKernel(E, cand));
        catch err
            continue;
        end try;
    end for;
    return #codomains gt 0, codomains;
end intrinsic;

//////////////////////////////////////////////////////////////////////////
// Accessors (Denotational semantics / Accessor formulas)
//////////////////////////////////////////////////////////////////////////

intrinsic MayBeReducible(ell::RngIntElt, primes::SeqEnum, info::Rec) -> BoolElt
{Whether ell may lie in R(E), given the sequence primes and IsogenyPrimesInfo
record info returned by IsogenyPrimes(E).}
    require info`Source eq "IsogenyPrimes" :
        "MayBeReducible: info must be an IsogenyPrimesInfo record (Source mismatch)";
    require IsPrime(ell) : "MayBeReducible: ell must be prime";
    if info`Kind eq "CMFamily" then
        famClause := (info`CMConductor mod ell ne 0)
            and (KroneckerSymbol(info`CMFundamentalDiscriminant, ell) ne -1);
        return famClause or (ell in primes);
    end if;
    // Kind eq "Finite"
    return ell in primes;
end intrinsic;

intrinsic MayBeCongruent(ell::RngIntElt, primes::SeqEnum, info::Rec) -> BoolElt
{Whether ell may be a congruence prime between E1 and E2, given the sequence
primes and CongruencePrimesInfo record info returned by CongruencePrimes(E1, E2).}
    require info`Source eq "CongruencePrimes" :
        "MayBeCongruent: info must be a CongruencePrimesInfo record (Source mismatch)";
    require IsPrime(ell) : "MayBeCongruent: ell must be prime";
    if info`Kind eq "AllPrimes" or info`Kind eq "Undecided" then
        return true;
    end if;
    // Kind eq "Finite"
    return ell in primes;
end intrinsic;

//////////////////////////////////////////////////////////////////////////
// Billerey reducible-prime machinery (branch 2: absolute degree >= 2, no
// geometric CM).
//
// Port of Sage 10.8 gal_reps_number_field.py (Billerey_P_l / _B_l / _R_q /
// _B_bound / _R_bound / Frobenius_filter / reducible_primes_Billerey), with
// the design spec's pinned deviations recorded inline. Star/power operations
// use CHIMP:
//   PowerCharacteristicPolynomial(f, n)  roots -> n-th powers  (adams_operator)
//   TensorCharacteristicPolynomial(f, g) roots -> products a*b (composed_op mul)
//////////////////////////////////////////////////////////////////////////

// Fetch a CHIMP intrinsic by name at call time. Nothing in this package's call
// graph may name a CHIMP intrinsic statically: on Magma 2.29-8 an absent one
// aborts at bind time, before the top-of-body requires can report it. This
// error is the backstop for callers that slip past those requires.
function ChimpIntrinsic(name)
    ok, f := IsIntrinsic(name);
    error if not ok, Sprintf(
        "IsogenyPrimes package: CHIMP is not attached (%o absent); AttachSpec CHIMP first", name);
    return f;
end function;

// Composed product (Sage composed_op with mul) of a nonempty polynomial seq.
function StarProductPolys(fs)
    tcp := ChimpIntrinsic("TensorCharacteristicPolynomial");
    g := fs[1];
    for i in [2 .. #fs] do
        g := tcp(g, fs[i]);
    end for;
    return g;
end function;

// Sage compose_power(f, k): degree-(deg f)^k polynomial whose roots are all
// k-fold products of roots of f. k = 0 is (x - 1); k >= 1 is the star of k
// copies of f.
function ComposePower(f, k)
    if k eq 0 then
        return Parent(f).1 - 1;
    end if;
    return StarProductPolys([ f : i in [1 .. k] ]);
end function;

// Residue characteristic of a conductor-support entry (the rational prime
// below a prime ideal; the integer itself over Q).
function ResidueCharacteristic(P)
    if Type(P) eq RngIntElt then
        return P;
    end if;
    return PrimeDivisors(Norm(P))[1];
end function;

// Auxiliary-prime admissibility (BOTH phases): r >= 5 prime not dividing
// 6 * Disc(K) * Norm(cond E). The conductor-based rule (not Sage's model
// discriminant) is the spec's second pinned deviation. r | Norm(cond E) iff r
// is a residue characteristic of the conductor support, precomputed in badRat;
// r | 6 iff r < 5 for a prime.
function IsAdmissiblePrime(r, DiscK, badRat)
    return r ge 5 and IsPrime(r) and (DiscK mod r ne 0) and (r notin badRat);
end function;

// Billerey P_l^* (eq (9)): star over primes q | l of
// PowerCharacteristicPolynomial(FrobeniusCharpoly(E, q), 12 e_q). Returns
// <false, _> when the model (globally integral, minimal where possible) is
// non-minimal at some q | l so FrobeniusCharpoly would reject it; the phase
// then skips l (safe: fewer auxiliary primes still give a superset). On a model
// minimal at good primes -- which admissibility guarantees -- this never fires.
function BilPlStar(E, l, DiscE)
    pcp := ChimpIntrinsic("PowerCharacteristicPolynomial");
    OK := Integers(BaseRing(E));
    factors := [];
    for pr in Decomposition(OK, l) do
        q := pr[1];
        if Valuation(DiscE, q) ne 0 then
            return false, _;
        end if;
        cp := FrobeniusCharpoly(E, q);
        Append(~factors, pcp(cp, 12 * pr[2]));
    end for;
    return true, StarProductPolys(factors);
end function;

// Sage Billerey_B_l with accumulator B (B = 0 gives the true B_l): product over
// k = 0 .. d/2 of GCD(P_l^*(l^{12k}), B), or 0 as soon as a factor vanishes.
// GCD(x, 0) = |x|, so values are nonnegative (matching Sage's ZZ.gcd) and the
// inert gate compares like-signed quantities. Second return is false when l was
// skipped (non-minimal model at some q | l).
function BilBlValue(E, l, d, B, DiscE)
    ok, P := BilPlStar(E, l, DiscE);
    if not ok then
        return 0, false;
    end if;
    Bl := 1;
    for k in [0 .. d div 2] do
        factor := Evaluate(P, l^(12*k));
        if factor eq 0 then
            return 0, true;
        end if;
        Bl *:= GCD(factor, B);
    end for;
    return Bl, true;
end function;

// Sage Billerey_R_q with accumulator B and the pinned deviation h = ord([q]):
// P = Power(charpoly_q, 12 h); Q = Power(minpoly(gamma), 12) with gamma a
// generator of the principal q^h; product over k = 0 .. d/2 of
// GCD(Res(P, compose_power(Q, k)), B), the k = 0 factor included (the term the
// PQM file drops), or 0 as soon as a resultant vanishes.
function BilRqValue(E, q, d, h, gamma, B)
    pcp := ChimpIntrinsic("PowerCharacteristicPolynomial");
    Zx := PolynomialRing(Integers());
    P := pcp(FrobeniusCharpoly(E, q), 12*h);
    Qpoly := pcp(Zx ! MinimalPolynomial(gamma), 12);
    Rq := 1;
    for k in [0 .. d div 2] do
        factor := Resultant(P, ComposePower(Qpoly, k));
        if factor eq 0 then
            return 0;
        end if;
        Rq *:= GCD(factor, B);
    end for;
    return Rq;
end function;

intrinsic BillereyBl(E::CrvEll, l::RngIntElt) -> RngIntElt
{Billerey's B_l (Theorem 2.4 of arXiv:0908.1084) for E over a number field at
the rational prime l: the product over k = 0 .. d/2 of P_l^* evaluated at
l^(12k), where P_l^* is equation (9). Exposed for the inert-principal gate.}
    require IsIntrinsic("PowerCharacteristicPolynomial") :
        "BillereyBl: CHIMP is not attached (PowerCharacteristicPolynomial absent); AttachSpec CHIMP first";
    require Type(BaseRing(E)) eq FldNum :
        "BillereyBl: E must be defined over a number field";
    d := AbsoluteDegree(BaseRing(E));
    require d ge 2 : "BillereyBl: Billerey's theorem requires absolute degree >= 2";
    require IsPrime(l) : "BillereyBl: l must be a rational prime";
    b := BilBlValue(E, l, d, 0, Discriminant(E));
    return b;
end intrinsic;

intrinsic BillereyRq(E::CrvEll, q::RngOrdIdl) -> RngIntElt
{Billerey's R_q (Theorem 2.8 of arXiv:0908.1084) for E over a number field at
the prime ideal q, with the pinned deviation h = ord([q]) in Cl(K) in place of
the class number. The full product over k = 0 .. d/2, k = 0 included. Exposed
for the inert-principal gate (R_q = B_l for l inert, q = (l)).}
    require IsIntrinsic("PowerCharacteristicPolynomial") :
        "BillereyRq: CHIMP is not attached (PowerCharacteristicPolynomial absent); AttachSpec CHIMP first";
    require Type(BaseRing(E)) eq FldNum :
        "BillereyRq: E must be defined over a number field";
    require IsPrime(q) : "BillereyRq: q must be a prime ideal";
    K := BaseRing(E);
    d := AbsoluteDegree(K);
    require d ge 2 : "BillereyRq: Billerey's theorem requires absolute degree >= 2";
    ClK, mClK := ClassGroup(K);
    h := Order(q @@ mClK);
    _, gamma := IsPrincipal(q^h);
    return BilRqValue(E, q, d, h, gamma, 0);
end intrinsic;

// Candidate-set prime factors of a running gcd, or the formal TOP sentinel
// while the accumulator is still 0 (never equal to any finite set, so the
// plateau rule never fires on it; PrimeDivisors is never called on 0).
function CandidateSet(G)
    if G eq 0 then
        return {Integers()|}, true;
    end if;
    return Seqset(PrimeDivisors(G)), false;
end function;

// B-phase with auto-escalation. Accumulates B := GCD(B, B_r) over admissible r
// up to a doubling bound (AuxBound .. MaxAuxBound). Returns
// <B, finalBound, stabilized, succeeded>; succeeded is false exactly when B is
// still 0 at the cap (the caller then runs the R-phase).
function RunBPhase(E, d, DiscE, DiscK, badRat, AuxBound, MaxAuxBound)
    B := 0;
    bound := AuxBound;
    prevSet := {Integers()|};
    prevIsTop := true;
    lastR := 4;
    while true do
        r := NextPrime(lastR);
        while r le bound do
            if IsAdmissiblePrime(r, DiscK, badRat) then
                bl, used := BilBlValue(E, r, d, B, DiscE);
                if used then
                    B := GCD(B, bl);
                    if B eq 1 then
                        return B, bound, true, true;
                    end if;
                end if;
            end if;
            lastR := r;
            r := NextPrime(r);
        end while;
        curSet, curIsTop := CandidateSet(B);
        vprintf IsogenyPrimes, 2: "  B-phase bound %o: support %o\n",
            bound, curIsTop select "TOP" else Sprint(SetToSequence(curSet));
        if bound eq MaxAuxBound then
            if AuxBound eq MaxAuxBound then
                stab := (B eq 1);
            else
                stab := (not curIsTop) and (not prevIsTop) and (curSet eq prevSet);
            end if;
            return B, bound, stab, (B ne 0);
        end if;
        if (not curIsTop) and (not prevIsTop) and (curSet eq prevSet) then
            return B, bound, true, true;
        end if;
        prevSet := curSet;
        prevIsTop := curIsTop;
        bound := Min(2*bound, MaxAuxBound);
    end while;
end function;

// R-phase with the same escalation over prime ideals q above admissible
// rationals, increasing norm (Decomposition order within a rational prime).
// h = ord([q]); gamma a generator of q^h. Each norm band (prevBound, bound] is
// visited once. Returns <R, finalBound, stabilized>; R = 0 at the cap means the
// R-phase also failed (the caller errors).
function RunRPhase(E, d, DiscE, DiscK, badRat, ClK, mClK, AuxBound, MaxAuxBound)
    OK := Integers(BaseRing(E));
    R := 0;
    bound := AuxBound;
    prevBound := 0;
    prevSet := {Integers()|};
    prevIsTop := true;
    while true do
        for r in PrimesInInterval(5, bound) do
            if not IsAdmissiblePrime(r, DiscK, badRat) then
                continue;
            end if;
            for pr in Decomposition(OK, r) do
                q := pr[1];
                nq := Norm(q);
                if nq gt bound or nq le prevBound then
                    continue;
                end if;
                if Valuation(DiscE, q) ne 0 then
                    continue;
                end if;
                h := Order(q @@ mClK);
                _, gamma := IsPrincipal(q^h);
                rq := BilRqValue(E, q, d, h, gamma, R);
                if rq ne 0 then
                    R := GCD(R, rq);
                    if R eq 1 then
                        return R, bound, true;
                    end if;
                end if;
            end for;
        end for;
        curSet, curIsTop := CandidateSet(R);
        vprintf IsogenyPrimes, 2: "  R-phase bound %o: support %o\n",
            bound, curIsTop select "TOP" else Sprint(SetToSequence(curSet));
        if bound eq MaxAuxBound then
            if AuxBound eq MaxAuxBound then
                stab := (R eq 1);
            else
                stab := (not curIsTop) and (not prevIsTop) and (curSet eq prevSet);
            end if;
            return R, bound, stab;
        end if;
        if (not curIsTop) and (not prevIsTop) and (curSet eq prevSet) then
            return R, bound, true;
        end if;
        prevSet := curSet;
        prevIsTop := curIsTop;
        prevBound := bound;
        bound := Min(2*bound, MaxAuxBound);
    end while;
end function;

// Frobenius filter (Sage Frobenius_filter, spec-pinned by FilterBound not
// patience): discard ell that has a good ideal q, Norm(q) <= FilterBound, with
// Frobenius charpoly x^2 - a x + N irreducible mod ell and constant term N
// nonzero mod ell. For odd ell that is KroneckerSymbol(a^2 - 4N, ell) = -1
// (which already forces ell not dividing N); ell = 2 is the a,N both-odd case.
// Genuine reducible primes never get a witness, so containment is preserved.
function FrobeniusFilterPrimes(E, cand, DiscE, FilterBound)
    OK := Integers(BaseRing(E));
    survivors := cand;
    if IsEmpty(survivors) then
        return survivors;
    end if;
    for r in PrimesUpTo(FilterBound) do
        for pr in Decomposition(OK, r) do
            q := pr[1];
            if Norm(q) gt FilterBound or Valuation(DiscE, q) ne 0 then
                continue;
            end if;
            cp := FrobeniusCharpoly(E, q);
            Nq := Integers() ! Coefficient(cp, 0);
            a := -(Integers() ! Coefficient(cp, 1));
            disc := a^2 - 4*Nq;
            toRemove := {Integers()|};
            for ell in survivors do
                if ell eq 2 then
                    if IsOdd(a) and IsOdd(Nq) then
                        Include(~toRemove, 2);
                    end if;
                elif Nq mod ell ne 0 and KroneckerSymbol(disc, ell) eq -1 then
                    Include(~toRemove, ell);
                end if;
            end for;
            survivors diff:= toRemove;
            if IsEmpty(survivors) then
                return survivors;
            end if;
        end for;
    end for;
    return survivors;
end function;

//////////////////////////////////////////////////////////////////////////
// CM branch machinery (branch 3: absolute degree >= 2, geometric CM).
//
// Single-implementation (the clean-room stream is exempt from CM), protected
// by the construction oracle and fixtures. Faithful port of Sage 10.8
// isogeny_class.py isogeny_degrees_cm(E): a finite set of primes generating
// the CM isogeny class, obtained by class-group enumeration over the orders
// between End(E) = O (discriminant d) and the maximal order, then trimmed by
// the Frobenius filter. C_CM := that ported output; the design spec's branch-3
// construction folds the finite fuzz into L := Sort(C_CM join {p : p divides
// f * Disc(K) * Norm(cond E)}).
//////////////////////////////////////////////////////////////////////////

// Sage BinaryQF.small_prime_value: the smallest prime represented by the
// positive-definite form a x^2 + b xy + c y^2, substituting x in [-B, B) and
// y in [0, B) with B growing by 10 (an elementary value search, matching the
// reference verbatim). Only ever used for the horizontal (same-order) primes,
// which arise solely under rational CM and are then split primes already in
// F(D_F, f); the specific value therefore never affects the CMFamily
// denotation (F covers them), so exact agreement with PARI's generator choice
// is not required here.
function SmallPrimeValueOfForm(a, b, c)
    B := 10;
    while true do
        vals := {Integers()|};
        for x in [-B .. B - 1] do
            for y in [0 .. B - 1] do
                v := a*x^2 + b*x*y + c*y^2;
                if v gt 1 and IsPrime(v) then
                    Include(~vals, v);
                end if;
            end for;
        end for;
        if #vals gt 0 then
            return Min(vals);
        end if;
        error if B ge 1000, "SmallPrimeValueOfForm: no prime value found below the bound";
        B +:= 10;
    end while;
end function;

// Port of Sage gal_reps_number_field.py Frobenius_filter (patience = 100). E
// MUST already be a global integral model (the CM branch passes the normalized
// model). 2 is handled first: removed, then re-added iff the 2-division
// polynomial is reducible (a K-rational 2-torsion x-coordinate, i.e. a
// possible 2-isogeny). For odd ell, ell is discarded as soon as some good
// prime P -- taken in rational-prime order over the first 100 good ideals --
// has Frobenius discriminant a^2 - 4 N(P) a nonsquare mod ell (Legendre = -1),
// which certifies E[ell] irreducible there. A genuine reducible prime never
// acquires such a witness, so the filter only removes provable non-isogeny
// primes; matching Sage's patience bound keeps C_CM equal to its output.
function SageFrobeniusFilter(E, Lset)
    K := BaseRing(E);
    N := Conductor(E);
    L := [ l : l in Sort(SetToSequence(Lset)) | l ne 2 ];
    include2 := false;
    if 2 in Lset then
        include2 := not IsIrreducible(DivisionPolynomial(E, 2));
    end if;
    numP := 0;
    p := 2;
    while numP lt 100 and #L gt 0 do
        if Type(K) eq FldRat then
            ideals := [ Integers() | p ];
        else
            ideals := [ pr[1] : pr in Decomposition(Integers(K), p) ];
        end if;
        for P in ideals do
            if numP ge 100 then
                break;
            end if;
            if Valuation(N, P) ne 0 then
                continue;               // bad reduction: Sage's iterator skips it
            end if;
            if Type(K) eq FldRat then
                a := TraceOfFrobenius(ChangeRing(E, GF(p)));
                Nq := p;
            else
                // E need not be locally minimal at P even though P is good
                // (conductor-based skip above): over class number > 1 fields
                // NormalizedModel can only return a globally integral model,
                // and Reduction requires local minimality at P. MinimalModel(E,P)
                // always succeeds for a single prime, regardless of class number.
                a := TraceOfFrobenius(Reduction(MinimalModel(E, P), P));
                Nq := Norm(P);
            end if;
            disc := a^2 - 4*Nq;
            L := [ l : l in L | KroneckerSymbol(disc, l) ne -1 ];
            numP +:= 1;
            if #L eq 0 then
                break;
            end if;
        end for;
        p := NextPrime(p);
    end while;
    result := Seqset(L);
    if include2 then
        Include(~result, 2);
    end if;
    return result;
end function;

// Port of isogeny_degrees_cm(E) returning the prime set C_CM (post filter). d
// is the CM ORDER discriminant (= f^2 * D_F from HasComplexMultiplication);
// rationalCM = E.has_rational_cm() (the CM order is defined over K, tested by
// the caller as IsSquare(K ! D_F)). Correspondence with the Sage reference:
//   n            base_field().absolute_degree(), doubled if not rationalCM,
//                times an extra unit factor for d = -3 (x3) and d = -4 (x2);
//   Divisors(n)  Sage's divs = n.divisors();
//   ClassNumber  data[0] from pari(d).quadclassunit() -> n/(2h);
//   ClassGroup(QuadraticForms(d)) generators + SmallPrimeValueOfForm
//                Qs = [BinaryQF(q) for q in data[2]] + Q.small_prime_value()
//                (horizontal same-order primes, rational CM only).
// Vertical primes: odd ramified (all if not rationalCM; else the upward l^2|d
// and the downward l | n/2h), downward split (kronecker +1), downward inert
// (kronecker -1). Initial {2} (plus 3 for d = -3).
function CMIsogenyDegrees(E, d, rationalCM)
    K := BaseRing(E);
    n := AbsoluteDegree(K);
    if not rationalCM then
        n *:= 2;
    end if;
    if d eq -4 then
        n *:= 2;
    end if;
    if d eq -3 then
        n *:= 3;
    end if;
    divs := Divisors(n);

    L := (d eq -3) select {Integers() | 2, 3} else {Integers() | 2};

    // Odd prime divisors of d (Sage d.odd_part().prime_factors()).
    ram_l := [ p : p in PrimeDivisors(Abs(d)) | IsOdd(p) ];

    // The class group of the order of discriminant d (imaginary quadratic,
    // possibly non-maximal); needed for n/2h and, if rationalCM, the horizontal
    // generators. Computed once here.
    if rationalCM then
        Forms := QuadraticForms(d);
        A, mA := ClassGroup(Forms);
        h := #A;
        n_over_2h := n div (2*h);
    end if;

    // (a) Ramified primes.
    if not rationalCM then
        L join:= Seqset(ram_l);
    else
        // Upward primes (index divided by l): l^2 | d.
        L join:= { l : l in ram_l | Valuation(Abs(d), l) gt 1 };
        // Downward ramified primes: l | n/2h.
        L join:= { l : l in ram_l | n_over_2h mod l eq 0 };
    end if;

    // (b) Downward split primes: l = lm1 + 1 prime, split (kronecker +1).
    L join:= { lm1 + 1 : lm1 in divs
               | IsPrime(lm1 + 1) and KroneckerSymbol(d, lm1 + 1) eq 1 };

    // (c) Downward inert primes: l = lp1 - 1 prime, inert (kronecker -1).
    L join:= { lp1 - 1 : lp1 in divs
               | IsPrime(lp1 - 1) and KroneckerSymbol(d, lp1 - 1) eq -1 };

    // Horizontal primes (rational CM only): one prime per class-group generator.
    if rationalCM then
        for g in Generators(A) do
            a, b, c := Explode(Eltseq(mA(g)));
            Include(~L, SmallPrimeValueOfForm(a, b, c));
        end for;
    end if;

    return SageFrobeniusFilter(E, L);
end function;

//////////////////////////////////////////////////////////////////////////
// Main intrinsics
//
// Task 4 scaffolds signatures, parameter preconditions, and model
// normalization only; the branch algorithms land in later tasks.
//////////////////////////////////////////////////////////////////////////

intrinsic IsogenyPrimes(E::CrvEll :
    AuxBound := 200, MaxAuxBound := 1600,
    FilterBound := 1000) -> SeqEnum, Rec
{A certified superset L of R(E), the primes ell such that E[ell] is
reducible as a Gal(Qbar/K)-module, plus a tagged IsogenyPrimesInfo record.
See the design spec's Denotational semantics section for the exact
guarantee.}
    require Type(BaseRing(E)) in {FldRat, FldNum} :
        "IsogenyPrimes: E must be defined over the rationals or a number field";
    require 0 lt AuxBound and AuxBound le MaxAuxBound :
        "IsogenyPrimes: require 0 < AuxBound <= MaxAuxBound";
    require FilterBound gt 0 : "IsogenyPrimes: require FilterBound > 0";

    Enorm := NormalizedModel(E);
    bad := ConductorSupport(Enorm);
    vprint IsogenyPrimes, 1: "IsogenyPrimes: normalized model", Enorm, "; bad primes", bad;

    K := BaseRing(Enorm);
    if AbsoluteDegree(K) eq 1 then
        // Branch 1 (checked first, CM or not): Kind := "Finite", Exact := true.
        // L is exactly Mazur's list filtered by HasPrimeIsogeny; no class-level
        // pre-pass and bound parameters are ignored (they are meaningless here).
        MAZUR := [2, 3, 5, 7, 11, 13, 17, 19, 37, 43, 67, 163];
        L := [ p : p in MAZUR | HasPrimeIsogeny(Enorm, p) ];
        info := rec< IsogenyPrimesInfo |
            Source := "IsogenyPrimes",
            Kind := "Finite",
            Exact := true,
            Stabilized := true,
            BoundsUsed := [Integers()|] >;

        // IsCM / CM fields are reported regardless of branch (record contract).
        isCM, D := HasComplexMultiplication(Enorm);
        info`IsCM := isCM;
        if isCM then
            Df := FundamentalDiscriminant(D);
            ok, f := IsSquare(D div Df);
            error if not ok, "IsogenyPrimes: internal error deriving the CM conductor from the order discriminant";
            info`CMOrderDiscriminant := D;
            info`CMFundamentalDiscriminant := Df;
            info`CMConductor := f;
            // CMInBaseField: pinned as IsSquare(K ! D_F); over Q (or any
            // degree-one field) this is always false since D_F < 0.
            info`CMInBaseField := IsSquare(K ! Df);
        end if;

        vprint IsogenyPrimes, 1: "IsogenyPrimes: branch 1 (absolute degree 1); L =", L;
        return L, info;
    end if;

    // --- Branches 2/3: absolute degree >= 2. Report CM status, dispatch on it. ---
    // HasComplexMultiplication is analytic; it returns geometric CM together
    // with the ORDER discriminant D_O = f^2 * D_F, and errors loudly if it
    // cannot decide (we never guess CM status).
    isCM, D := HasComplexMultiplication(Enorm);
    info := rec< IsogenyPrimesInfo |
        Source := "IsogenyPrimes",
        Kind := "Finite",
        Exact := false,
        IsCM := isCM >;
    if isCM then
        // --- Branch 3: degree >= 2 geometric CM. Kind := "CMFamily" iff the CM
        // is defined over K (CMInBaseField), else "Finite"; Exact stays false;
        // Stabilized := true; BoundsUsed := []. L := C_CM (ported
        // isogeny_degrees_cm) union the fuzz f * Disc(K) * Norm(cond E). ---
        Df := FundamentalDiscriminant(D);
        ok, f := IsSquare(D div Df);
        error if not ok,
            "IsogenyPrimes: internal error deriving the CM conductor from the order discriminant";
        // CMInBaseField: all geometric endomorphisms are defined over K iff the
        // CM field Q(sqrt(D_F)) embeds in K, tested exactly as IsSquare(K ! D_F).
        rationalCM := IsSquare(K ! Df);
        info`CMOrderDiscriminant := D;
        info`CMFundamentalDiscriminant := Df;
        info`CMConductor := f;
        info`CMInBaseField := rationalCM;
        info`Kind := rationalCM select "CMFamily" else "Finite";
        info`Stabilized := true;
        info`BoundsUsed := [Integers() | ];

        C_CM := CMIsogenyDegrees(Enorm, D, rationalCM);
        DiscK := Discriminant(Integers(K));
        fuzz := Seqset(PrimeDivisors(Abs(f * DiscK * Norm(Conductor(Enorm)))));
        L := Sort(SetToSequence(C_CM join fuzz));
        vprint IsogenyPrimes, 1:
            "IsogenyPrimes: branch 3 (geometric CM, Kind", info`Kind, "); C_CM",
            Sort(SetToSequence(C_CM)), "fuzz", Sort(SetToSequence(fuzz)), "-> L", L;
        return L, info;
    end if;

    // --- Branch 2: degree >= 2, no geometric CM (certified superset of R(E)) ---
    require IsIntrinsic("PowerCharacteristicPolynomial") :
        "IsogenyPrimes: CHIMP is not attached (PowerCharacteristicPolynomial absent); AttachSpec CHIMP first";

    OK := Integers(K);
    DiscE := Discriminant(Enorm);
    DiscK := Discriminant(OK);
    d := AbsoluteDegree(K);
    badRat := { ResidueCharacteristic(P) : P in bad };
    vprint IsogenyPrimes, 1:
        "IsogenyPrimes: branch 2 (abs degree", d, "non-CM); bad rational primes", badRat;

    B, bBound, stabilized, bOK := RunBPhase(Enorm, d, DiscE, DiscK, badRat, AuxBound, MaxAuxBound);
    if bOK then
        G := B;
        finalRBound := 0;
        vprint IsogenyPrimes, 1: "IsogenyPrimes: B-phase support", PrimeDivisors(G), "at bound", bBound;
    else
        vprint IsogenyPrimes, 1: "IsogenyPrimes: B-phase vanished; running R-phase";
        ClK, mClK := ClassGroup(K);
        R, rBound, stabilized := RunRPhase(Enorm, d, DiscE, DiscK, badRat, ClK, mClK, AuxBound, MaxAuxBound);
        error if R eq 0,
            "IsogenyPrimes: both Billerey phases vanished at MaxAuxBound; geometric CM was likely missed or MaxAuxBound is too small";
        G := R;
        bBound := MaxAuxBound;
        finalRBound := rBound;
        vprint IsogenyPrimes, 1: "IsogenyPrimes: R-phase support", PrimeDivisors(G), "at bound", rBound;
    end if;

    // Conservative candidate enlargement, then the Frobenius filter.
    cand := Seqset(PrimeDivisors(G)) join {2, 3}
            join Seqset(PrimeDivisors(Abs(DiscK))) join badRat;
    L := FrobeniusFilterPrimes(Enorm, cand, DiscE, FilterBound);
    Lsorted := Sort(SetToSequence(L));

    info`Stabilized := stabilized;
    info`BoundsUsed := [ bBound, finalRBound, FilterBound ];
    vprint IsogenyPrimes, 1:
        "IsogenyPrimes: candidates", Sort(SetToSequence(cand)), "-> filtered L", Lsorted;
    return Lsorted, info;
end intrinsic;

//////////////////////////////////////////////////////////////////////////
// CongruencePrimes machinery (flow steps 2 and 4).
//////////////////////////////////////////////////////////////////////////

// Step 2 gcd loop. Over prime ideals q good for BOTH curves (conductor support
// after model normalization; rational primes over Q), in increasing norm,
// accumulate G := GCD of p_q * (a_q(E1) - a_q(E2)), p_q the residue
// characteristic. The p_q factor keeps ell = p_q admissible at its own q, so
// the survivors are exactly the primes ell with a_q(E1) = a_q(E2) mod ell at
// every sampled q. Zero-sentinel (G = 0 is the formal TOP, never PrimeFactored,
// never plateau-equal to a finite set) and the doubling escalation are the
// B-phase's, including the single-evaluation Stabilized rule. Returns
// <G, finalBound, stabilized>; G = 0 at the cap means the traces never
// separated (the caller decides Undecided vs the BFS).
function RunCongruenceGcdPhase(E1, E2, bad1, bad2, NormBound, MaxNormBound)
    K := BaseRing(E1);
    isQ := Type(K) eq FldRat;
    if not isQ then
        OK := Integers(K);
    end if;
    bad1set := Set(bad1);
    bad2set := Set(bad2);
    G := 0;
    bound := NormBound;
    prevBound := 0;
    prevSet := {Integers()|};
    prevIsTop := true;
    while true do
        for p in PrimesUpTo(bound) do
            if isQ then
                qs := [ Integers() | p ];
            else
                qs := [ pr[1] : pr in Decomposition(OK, p) ];
            end if;
            for q in qs do
                nq := isQ select p else Norm(q);
                if nq gt bound or nq le prevBound then
                    continue;               // outside this norm band
                end if;
                if q in bad1set or q in bad2set then
                    continue;               // bad for one of the curves
                end if;
                a1 := - (Integers() ! Coefficient(FrobeniusCharpoly(E1, q), 1));
                a2 := - (Integers() ! Coefficient(FrobeniusCharpoly(E2, q), 1));
                val := p * (a1 - a2);
                if val ne 0 then
                    G := GCD(G, val);
                    if G eq 1 then
                        return G, bound, true;       // absorbing (empty L)
                    end if;
                end if;
            end for;
        end for;
        curSet, curIsTop := CandidateSet(G);
        vprintf IsogenyPrimes, 2: "  congruence gcd bound %o: support %o\n",
            bound, curIsTop select "TOP" else Sprint(SetToSequence(curSet));
        if bound eq MaxNormBound then
            if NormBound eq MaxNormBound then
                stab := (G eq 1);
            else
                stab := (not curIsTop) and (not prevIsTop) and (curSet eq prevSet);
            end if;
            return G, bound, stab;
        end if;
        if (not curIsTop) and (not prevIsTop) and (curSet eq prevSet) then
            return G, bound, true;
        end if;
        prevSet := curSet;
        prevIsTop := curIsTop;
        prevBound := bound;
        bound := Min(2*bound, MaxNormBound);
    end while;
end function;

// Step 4 (degree >= 2) certification: BFS from E1 through prime-degree isogenies
// (HasPrimeIsogeny validated kernels, degree <= primeBound, depth <= maxDepth),
// expanding over ALL validated codomains per prime degree -- a single-witness
// expansion could miss a reachable E2. A node certifies iff
// IsIsomorphic(node, E2) over K; equal j-invariants are only a prefilter (a
// nontrivial quadratic twist shares E2's j but twists the representation, so
// bare j-equality must never certify). Nodes are deduped by K-isomorphism
// within equal-j buckets (visited maps j -> pairwise non-isomorphic
// representatives). Returns whether some node is K-isomorphic to E2; the search
// can only confirm, never refute.
function CertifyByIsogenyBFS(E1, E2, primeBound, maxDepth)
    jTarget := jInvariant(E2);
    primes := PrimesUpTo(primeBound);
    visited := AssociativeArray();

    jE1 := jInvariant(E1);
    visited[jE1] := [ E1 ];
    if jE1 eq jTarget and IsIsomorphic(E1, E2) then
        return true;                        // depth-0 node already matches
    end if;

    frontier := [ E1 ];
    depth := 0;
    while depth lt maxDepth and #frontier gt 0 do
        next := [];
        for C in frontier do
            for p in primes do
                _, codomains := HasPrimeIsogeny(C, p);
                for D in codomains do
                    jD := jInvariant(D);
                    isNew := true;
                    if IsDefined(visited, jD) then
                        for rep in visited[jD] do
                            if IsIsomorphic(D, rep) then
                                isNew := false;    // K-isomorphic duplicate
                                break;
                            end if;
                        end for;
                        if isNew then
                            visited[jD] := Append(visited[jD], D);
                        end if;
                    else
                        visited[jD] := [ D ];
                    end if;
                    if isNew then
                        if jD eq jTarget and IsIsomorphic(D, E2) then
                            return true;
                        end if;
                        Append(~next, D);
                    end if;
                end for;
            end for;
        end for;
        frontier := next;
        depth +:= 1;
    end while;
    return false;
end function;

intrinsic CongruencePrimes(E1::CrvEll, E2::CrvEll :
    NormBound := 1000, MaxNormBound := 8000,
    KnownIsogenous := false,
    CertificationPrimeBound := 100, CertificationDepth := 3) -> SeqEnum, Rec
{A certified superset L of the primes ell such that E1[ell] and E2[ell] have
isomorphic semisimplifications over their common base field, plus a tagged
CongruencePrimesInfo record. See the design spec's CongruencePrimes
semantics section for the exact guarantee.}
    require Type(BaseRing(E1)) in {FldRat, FldNum} :
        "CongruencePrimes: E1 must be defined over the rationals or a number field";
    require Type(BaseRing(E2)) in {FldRat, FldNum} :
        "CongruencePrimes: E2 must be defined over the rationals or a number field";
    require 0 lt NormBound and NormBound le MaxNormBound :
        "CongruencePrimes: require 0 < NormBound <= MaxNormBound";
    require CertificationPrimeBound ge 2 :
        "CongruencePrimes: require CertificationPrimeBound >= 2";
    require CertificationDepth ge 0 :
        "CongruencePrimes: require CertificationDepth >= 0";

    E1norm := NormalizedModel(E1);
    E2norm := NormalizedModel(E2);
    require BaseRing(E1norm) eq BaseRing(E2norm) :
        "CongruencePrimes: E1 and E2 must be defined over a common base field";
    K := BaseRing(E1norm);
    deg := AbsoluteDegree(K);
    vprint IsogenyPrimes, 1: "CongruencePrimes: normalized models", E1norm, E2norm;

    // Step 0: a consumer-supplied isogeny relation certifies all primes.
    // Needed because IsIsogenous errors on CrvEll over FldNum on 2.29-8, so an
    // isogeny known to the caller (e.g. from an M2 handler) cannot always be
    // re-derived in-engine.
    if KnownIsogenous then
        vprint IsogenyPrimes, 1: "CongruencePrimes: KnownIsogenous supplied -> AllPrimes";
        return [Integers()|], rec< CongruencePrimesInfo |
            Source := "CongruencePrimes",
            Kind := "AllPrimes",
            Exact := true,
            Stabilized := true,
            BoundsUsed := [0, 0, 0],
            CertificationMethod := "Supplied" >;
    end if;

    // Step 1: over Q, IsIsogenous decides all primes at once. It is guarded to
    // FldRat because it errors on CrvEll over FldNum on this Magma version;
    // a degree-one FldNum therefore falls through to the gcd loop instead.
    if deg eq 1 and Type(K) eq FldRat and IsIsogenous(E1norm, E2norm) then
        vprint IsogenyPrimes, 1: "CongruencePrimes: IsIsogenous over Q -> AllPrimes";
        return [Integers()|], rec< CongruencePrimesInfo |
            Source := "CongruencePrimes",
            Kind := "AllPrimes",
            Exact := true,
            Stabilized := true,
            BoundsUsed := [0, 0, 0],
            CertificationMethod := "IsIsogenous" >;
    end if;

    // Step 2: gcd of p_q * (a_q(E1) - a_q(E2)) over shared-good q, escalating.
    bad1 := ConductorSupport(E1norm);
    bad2 := ConductorSupport(E2norm);
    vprint IsogenyPrimes, 1: "CongruencePrimes: bad primes", bad1, bad2;
    G, finalBound, stabilized :=
        RunCongruenceGcdPhase(E1norm, E2norm, bad1, bad2, NormBound, MaxNormBound);

    // Step 3: nonzero gcd -> finite candidate superset (uncertified).
    if G ne 0 then
        L := PrimeDivisors(G);
        vprint IsogenyPrimes, 1: "CongruencePrimes: G =", G, "-> L", L;
        return L, rec< CongruencePrimesInfo |
            Source := "CongruencePrimes",
            Kind := "Finite",
            Exact := false,
            Stabilized := stabilized,
            BoundsUsed := [finalBound, 0, 0],
            CertificationMethod := "None" >;
    end if;

    // Step 4: G = 0 at the cap (traces agreed at every sampled q).
    if deg eq 1 then
        // Reachable even after IsIsogenous returned false: a finite cap need
        // not contain a separating prime. The only safe upper bound is all
        // primes; report it as Undecided (callers short-circuit on Kind).
        vprint IsogenyPrimes, 1: "CongruencePrimes: G = 0 at cap, degree 1 -> Undecided";
        return [Integers()|], rec< CongruencePrimesInfo |
            Source := "CongruencePrimes",
            Kind := "Undecided",
            Exact := false,
            Stabilized := false,
            BoundsUsed := [MaxNormBound, 0, 0],
            CertificationMethod := "None" >;
    end if;

    // Degree >= 2: bounded explicit-isogeny certification.
    vprint IsogenyPrimes, 1:
        "CongruencePrimes: G = 0 at cap, degree", deg, "-> isogeny BFS";
    found := CertifyByIsogenyBFS(E1norm, E2norm, CertificationPrimeBound, CertificationDepth);
    if found then
        vprint IsogenyPrimes, 1: "CongruencePrimes: BFS certified an isogeny -> AllPrimes";
        return [Integers()|], rec< CongruencePrimesInfo |
            Source := "CongruencePrimes",
            Kind := "AllPrimes",
            Exact := true,
            Stabilized := true,
            BoundsUsed := [MaxNormBound, CertificationPrimeBound, CertificationDepth],
            CertificationMethod := "ExplicitIsogeny" >;
    end if;
    vprint IsogenyPrimes, 1: "CongruencePrimes: BFS did not reach E2 within bounds -> Undecided";
    return [Integers()|], rec< CongruencePrimesInfo |
        Source := "CongruencePrimes",
        Kind := "Undecided",
        Exact := false,
        Stabilized := false,
        BoundsUsed := [MaxNormBound, CertificationPrimeBound, CertificationDepth],
        CertificationMethod := "None" >;
end intrinsic;
