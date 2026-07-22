// Shared support for the isogeny-primes test drivers
// (tests/test_isogenyprimes.m and tests/run_differential.m): the
// attach/guard preamble, the parser for the committed data file
// tests/corpus.txt (schema documented in its header), and the field/curve
// constructors. Both drivers load this file first, from the repository root.

// spec:="" is the red-state syntactic check: no engine spec, CHIMP not
// needed; the engine-presence guard below then quits 1.
if assigned spec then useSpec := spec; else useSpec := "magma/spec"; end if;
if useSpec ne "" then
    // CHIMP supplies the star/bracket charpoly intrinsics the engine needs.
    ok, _ := IsIntrinsic("PowerCharacteristicPolynomial");
    if not ok then
        printf "SUITE FAILED: CHIMP is not attached (PowerCharacteristicPolynomial absent); AttachSpec your CHIMP.spec first\n";
        quit 1;
    end if;
    AttachSpec(useSpec);
end if;
if not assigned section then section := "all"; end if;
if not assigned cmscope then cmscope := "1"; end if;
// Engine-presence guard: without the engine intrinsics the driver procedures
// fail at BIND time, which try/catch cannot intercept, so quit with an honest
// verdict before any of that noise.
okEngine, _ := IsIntrinsic("IsogenyPrimes");
okEngine2, _ := IsIntrinsic("CongruencePrimes");
if not (okEngine and okEngine2) then
    printf "SUITE FAILED: engine intrinsics not attached (red state)\n";
    quit 1;
end if;

R<x> := PolynomialRing(Rationals());

function BuildField(coeffs)
    if coeffs eq [0, 1] then return Rationals(); end if;
    // NumberField(R!coeffs) collapses a degree-1 defining polynomial to
    // FldRat in this Magma version; force a genuine FldNum here so the
    // degree-one dispatch path is actually exercised.
    if #coeffs eq 2 then return RationalsAsNumberField(); end if;
    return NumberField(R ! coeffs);
end function;

function Elt(K, v)
    if Type(K) eq FldRat then return K ! v[1]; end if;
    return &+[ (Rationals() ! v[i]) * K.1^(i-1) : i in [1..#v] ];
end function;

function BuildCurve(K, ainvs)
    return EllipticCurve([ Elt(K, a) : a in ainvs ]);
end function;

function PrimeAbove(K, p, i)
    return Decomposition(Integers(K), p)[i][1];
end function;

function HasPrefix(s, p)
    return #s ge #p and s[1..#p] eq p;
end function;

function IntSet(L)
    return { Integers() | p : p in L };
end function;

// ---------------------------------------------------------------------------
// tests/corpus.txt parsing. One record per data line; unassigned record
// fields mean the corpus field was empty (oracle dropped / not applicable /
// not pinned), which the drivers distinguish from a computed-empty list.
// ---------------------------------------------------------------------------

CorpusRF := recformat<
    id : MonStgElt, stratum : MonStgElt, ispair : BoolElt, deg1 : BoolElt,
    field, ainvs, ainvs2, cm, oE, soft, golden, dropped,
    congkind : MonStgElt, congstab : BoolElt, congbound : RngIntElt, congprimes,
    lowkind : MonStgElt, lowstab : BoolElt, lowbound : RngIntElt, lowprimes,
    regkind : MonStgElt, regexact : BoolElt, regstab : BoolElt, regprimes,
    regcert : MonStgElt, regcmdisc : RngIntElt >;

// Split on sep at bracket depth 0, preserving empty fields (Magma's Split
// drops them, which would shift every column after an absent value).
function SplitTop(s, sep)
    parts := [];
    cur := "";
    depth := 0;
    for i in [1..#s] do
        c := s[i];
        if c eq "[" then
            depth +:= 1;
        elif c eq "]" then
            depth -:= 1;
        end if;
        if c eq sep and depth eq 0 then
            Append(~parts, cur);
            cur := "";
        else
            cur cat:= c;
        end if;
    end for;
    Append(~parts, cur);
    return parts;
end function;

// The comma-separated tokens of one bracketed list, as strings; used for the
// dropped-oracle labels, which contain colons and are not Magma literals.
function StrList(s)
    error if #s lt 2 or s[1] ne "[" or s[#s] ne "]",
        "corpus: expected a bracketed list, got " cat s;
    if #s eq 2 then return []; end if;
    return SplitTop(s[2..#s-1], ",");
end function;

function CorpusRecord(line)
    f := SplitTop(line, ":");
    error if #f notin {14, 20},
        Sprintf("corpus line for %o has %o fields, want 14 or 20", f[1], #f);
    r := rec< CorpusRF | id := f[1], stratum := f[2] >;
    // Every numeric field is a Magma sequence literal, so eval parses it
    // (rationals included); only the dropped labels need StrList.
    r`field := eval f[3];
    r`ainvs := eval f[4];
    if #f eq 14 then
        r`ispair := false;
        r`deg1 := f[5] eq "1";
        if f[6] ne "" then r`cm := eval f[6]; end if;
        if f[7] ne "" then r`oE := eval f[7]; end if;
        if f[8] ne "" then r`soft := eval f[8]; end if;
        if f[9] ne "" then r`golden := eval f[9]; end if;
        r`dropped := StrList(f[10]);
        if f[11] ne "" then
            r`regkind := f[11];
            r`regexact := f[12] eq "1";
            r`regprimes := eval f[13];
            r`regcmdisc := StringToInteger(f[14]);
        end if;
    else
        r`ispair := true;
        r`ainvs2 := eval f[5];
        r`deg1 := f[6] eq "1";
        if f[7] ne "" then
            r`congkind := f[7];
            r`congstab := f[8] eq "1";
            r`congbound := StringToInteger(f[9]);
            if f[10] ne "" then r`congprimes := eval f[10]; end if;
        end if;
        if f[11] ne "" then
            r`lowkind := f[11];
            r`lowstab := f[12] eq "1";
            r`lowbound := StringToInteger(f[13]);
            if f[14] ne "" then r`lowprimes := eval f[14]; end if;
        end if;
        r`dropped := StrList(f[15]);
        if f[16] ne "" then
            r`regkind := f[16];
            r`regexact := f[17] eq "1";
            r`regstab := f[18] eq "1";
            r`regprimes := eval f[19];
            r`regcert := f[20];
        end if;
    end if;
    // deg1 is redundant with the field encoding; a mismatch means a data
    // line was corrupted by hand.
    error if r`deg1 ne (#r`field eq 2),
        Sprintf("corpus %o: deg1 flag contradicts the field degree", r`id);
    return r;
end function;

function LoadCorpus(path)
    entries := [* *];
    for line in Split(Read(path), "\n") do
        if #line eq 0 or line[1] eq "#" then continue; end if;
        Append(~entries, CorpusRecord(line));
    end for;
    error if #entries eq 0, "corpus has no data lines: " cat path;
    return entries;
end function;

// The classification the drivers share: an entry is CM iff the Sage CM
// oracle said so ([1,...]); a dropped CM oracle classifies as non-CM, the
// same direction the oracle side takes.
function IsCMEntry(r)
    return assigned r`cm and r`cm[1] eq 1;
end function;

// A dropped per-ell construction (oE:ell=<ell>) or Frobenius screen
// (oE:frob) makes the oE oracle a lower bound for this entry, so equality
// asserts must downgrade to containment.
function HasOEDrop(r)
    return exists{ d : d in r`dropped | HasPrefix(d, "oE:") };
end function;

corpusLoaded := false;
try
    corpus := LoadCorpus("tests/corpus.txt");
    corpusLoaded := true;
catch e
    printf "SUITE FAILED: cannot load tests/corpus.txt: %o\n", e`Object;
end try;
if not corpusLoaded then
    quit 1;
end if;
