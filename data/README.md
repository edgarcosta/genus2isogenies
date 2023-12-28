Each row of the output corresponds to an isogeny class.

This is the description of the columns:
1. Conductor
2. LMFDB isogeny class label if known, otherwise `NULL`
3. Trace hash, as in Section 4.3 of https://arxiv.org/abs/1602.03715
4. Modular invariants for the known curves in the isogeny class
5. Curve equations corresponding to the invariants in item 4
6. The primes ell that we considered to test for isogenies
7. Modular invariants for all the curves in the isogeny class (this is a superset of item 4)
8. Curve equations corresponding to the invariants in item 7 (up to isomorphism, this is a superset of item 5)
9. A matrix representing the weighted isogeny graph, where the weights are the square roots of the isogeny degrees.
Note that not every isogeny needs to be irreducible, i.e., some of the isogenies of degree p^4 may be obtained by the composition of two isogenies of degree p^2.
For example, all the 7^4 isogenies represented in `[[0,7,7],[7,0,49],[7,49,0]]` are decomposable into two isogenies of degree 7^2.

The input is given by columns 1 through 5.
