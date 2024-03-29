SetQuitOnError(true);
if false then
  CHIMP();
end if;
SetQuitOnError(false);

/*
Depends on CHIMP https://github.com/edgarcosta/CHIMP
From Magma/utils.m:
- ReplaceCharacter
- IndexFibers
From endomorphisms package:
- HyperellipticCurveExtra
- PeriodMatrix
- GeometricHomomorphismRepresentation
*/

// example how to run it
// time parallel -j 64 --eta -a data/output_g2database_2e20.txt.jun25 --joblog data/richelot_count.log  magma -b input:={} verify.m > data/richelot_count2
hwlpolys := 0 eq System("which hwlpolys > /dev/null"); // soon to be available in a new version of smalljac



function TracesOfFrobeniusQuick(C, B0, B1: exclude:={})
  if not IsIntegral(C) then C:=IntegralModel(C); end if;
  D := Integers()!Discriminant(C);
  filename := Sprintf("temp_%o", Getpid());
  if hwlpolys then
    cstr := StripWhiteSpace(Sprint(Coefficients(HyperellipticPolynomials(SimplifiedModel(C)))));
    k := Ceiling(Log(2,B1));
    cmd := Sprintf("hwlpolys %o %o 2 0 -1 1 %o", cstr, k, filename);
    _ := Pipe(cmd, "");
  else
    f := HyperellipticPolynomials(SimplifiedModel(C));
    _<x> := Parent(f);
    cstr := StripWhiteSpace(Sprint(f));
    cmd := Sprintf("lpdata %o \"%o\" %o 3", filename, cstr, B1);
    _ := Pipe(cmd, "");
    filename cat:= "_lpdata.txt";
  end if;
  S := [[StringToInteger(c):c in Split(r,",")]: r in Split(Read(filename)) | "," in r];
  System("rm " cat filename);
  X := IndexFibers(S,func<r|r[1]>:Unique,Project:=func<r|r[2]>);
  minp := 16*Genus(C)^2;
  res := [ p gt minp and IsDefined(X,p) select X[p] else p+1-#ChangeRing(C,GF(p)):p in PrimesInInterval(B0,B1) | not (p in exclude or IsDivisibleBy(D,p))];
  if not hwlpolys then
    res := [-elt : elt in res];
  end if;
  return res;
end function;

function RichelotClass(Cs)
    queue := [<elt, G2Invariants(elt)>: elt in Cs];
    res := [];
    while #queue gt 0 do
        //print [G2Invariants(elt) : elt in queue];
        Cp, g2p := Explode(queue[1]);
        queue := queue[2..#queue];
        if not &or[IsIsomorphic(elt[1], Cp) : elt in res | elt[2] eq g2p ] then
            Append(~res, <Cp, g2p>);
        else
            continue; // already took care of it
        end if;
        //if g2 in done then continue end if;
        for Cpp in RichelotIsogenousSurfaces(Cp) do
            g2pp := G2Invariants(Cpp);
            Append(~queue, <Cpp, g2pp>);
        end for;
    end while;
    return [IntegralModel(elt[1]) : elt in res];
end function;

if assigned debug then
  SetDebugOnError(true);
end if;

print_richelot := assigned richelot_count;
check_heuristic := not assigned noheuristic;
check_traces := not assigned notraces;

function Verify(input)
  s := Split(input, ":");
  Lhash := eval s[3];
  eqns := eval s[8];
  R<x> := PolynomialRing(Integers());
  M := eval s[9];
  M := Matrix([[elt^2 : elt in row] : row in M]);

  curves := [HyperellipticCurve(R!elt[1], R!elt[2]) : elt in eqns];

  richelotclass := RichelotClass(curves);
  if print_richelot then
    old_eqns := eval s[5];
    old_curves := [HyperellipticCurve(R!elt[1], R!elt[2]) : elt in old_eqns];
    print #{G2Invariants(elt) : elt in RichelotClass(old_curves)};
  end if;
  // Check that we obtained RichelotIsogenousSurfaces
  assert {G2Invariants(elt) : elt in richelotclass} subset {G2Invariants(c) : c in curves};

  if check_traces then
    // Check that aps match for p < 2^16, with some exceptions
    exclude := &join[SequenceToSet(PrimeDivisors(Integers()!Discriminant(elt))) : elt in curves];
    traces := {TracesOfFrobeniusQuick(elt, 1, 2^16 : exclude:=exclude) : elt in curves};
    assert #traces eq 1;
  end if;

  if check_heuristic then
    done := false;
    prec := 200;
    heuristicM := ZeroMatrix(Integers(), Nrows(M));
    while not done do
      curves := [HyperellipticCurveExtra(R!elt[1], R!elt[2], prec) : elt in eqns];
      periods := [PeriodMatrix(elt) : elt in curves];
      try
        for i->pi in periods do
          for j->pj in periods[1..i-1] do
            if M[i][j] eq 0 then continue; end if;
            // we only want rational
            heuristicM[i][j] := Integers()!Determinant(GeometricHomomorphismRepresentation(periods[i], periods[j], RationalsExtra(prec) : UpperBound:=1)[1][2]);
            heuristicM[j][i] := heuristicM[i][j];
          end for;
        end for;
        done := true;
      catch e
        prec := prec*2;
        heuristicM := ZeroMatrix(Integers(), Nrows(M));
      end try;
    end while;
    assert heuristicM eq M;
  end if;
  return true;
end function;

res := [Verify(ReplaceCharacter(elt, "\\", "")) : elt in Split(input, " ")];
exit;
