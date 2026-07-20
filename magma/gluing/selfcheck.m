intrinsic GluingSelfCheck() -> BoolElt
{Package smoke test: shims and basic invariants. Returns true or raises an assertion.}
    CC := ComplexField(60);
    tau := Matrix(CC, 2, 2, [CC.1, 0, 0, 2*CC.1]);  // already reduced, diagonal
    taured, g := SiegelReduce(tau);
    assert Nrows(taured) eq 2;
    ok, q := RecognizeRational(CC!(22/7) + CC!10^-50*CC.1);
    assert ok and q eq 22/7;
    ok, _ := RecognizeRational(CC!Pi(CC));
    assert not ok;
    return true;
end intrinsic;
