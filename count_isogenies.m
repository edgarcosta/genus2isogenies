// alternative script to count isogenies
/* Expected output:
2^2 419157
2^4 693519
3^2 11568
3^4 29742
5^2 415
5^4 2440
7^2 154
7^4 246
11^4 9
13^2 20
13^4 9
17^2 4
31^2 1
*/

data := [atoiii(Split(line, ":")[9]) :  line in Split(Read("output_g2database_2e20.txt"))];


// returns a Matrix M
function FilterMatrixOn(M, check)
    n := Nrows(M);
    m := Ncols(M);
    return Matrix(BaseRing(M), [[check(i,j) select M[i,j] else 0 : j in [1..m]] : i in [1..n]]);
end function;

function MinimalGraph(M)
    assert Transpose(M) eq M;
    // only prime entries
    is_prime := func<i,j | IsPrime(M[i,j])>;
    Mprime := FilterMatrixOn(M, is_prime);
    // everything reached in 2 steps via prime edges
    M2prime := Mprime^2;
    is_not_twostep := func<i,j | Mprime[i,j] ne 0 or M2prime[i,j] eq 0>;
    // keep only the edges that are either one-step or not reachable via combination of two one-steps
    return FilterMatrixOn(M, is_not_twostep);
end function;

counter := AssociativeArray();
for row in data do
    M := MinimalGraph2(Matrix(Integers(), row));
    for e->ct in Multiset(Eltseq(M)) do
        c := ct div 2;
        if IsDefined(counter, e) then
            counter[e] +:= c;
        else
            counter[e] := c;
        end if;
    end for;
end for;

for d in Sort(SetToSequence(Set([PrimeDivisors(elt)[1] : elt in Keys(counter) | elt ne 0]))) do
    if IsPrime(d) then
        if IsDefined(counter, d) then
            printf "%o^2 %o\n", d, counter[d];
        end if;
        if IsDefined(counter, d^2) then
            printf "%o^4 %o\n", d, counter[d^2];
        end if;
    end if;
end for;
