Each row of the output corresponds to an isogeny class.

This is the description of the columns:
1. Conductor
2. LMFDB isogeny class label if known, otherwise `NULL`
3. Trace hash, as in Section 4.3 of https://arxiv.org/abs/1602.03715
4. Modular invariants for the known curves in the isogeny class
5. Curve equations corresponding to the invariants in item 4
6. The primes ell that we considered to test for isogenies
7. Modular invariants for all the curves in the isogeny class (this is a supset of item 4)
8. Curve equations corresponding to the invariants in item 7 (up to isomormphism, this is a supset of item 5)
9. A matrix representing the weighted isogeny graph, where the weights are the square roots of the isogeny degrees

The input is given by columns 1 through 5.
