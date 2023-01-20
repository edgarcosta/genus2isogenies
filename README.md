# genus2isogenies


# Computation of isogeny classes
The computation of isogeny classes depends on SageMath and [pyhdme](https://github.com/edgarcosta/pyhdme)

Here are some examples:
```
$ sage -python genus2isogenies.py "704:704.a:2159575534599616393:[[-1856,129536,-44,-2948]]:[[[-1,-1,-2,1,1,-2],[1,1,1,1]]]"
704:704.a:2159575534599616393:[[-1856,129536,-44,-2948]]:[[[-1,-1,-2,1,1,-2],[1,1,1,1]]]:[2,3]:[[-1856,129536,-44,-2948],[771904,687149056,260876,-197548868]]:[[[0,4,4,-1,-2],[0,0,0,1]],[[-1140,-5166,-4646,2135,-236,4],[1]]]:[[0,9],[9,0]]
```

A different example using 12 threads (it searches for isogenies of degree 29^2 and 29^4):
```
sage -python genus2isogenies.py "976:976.a:78262366554954104:[[1012,52318,-976,-37088]]:[[[1,2,0,-2,-1,-1],[0,0,1,1]]]" 12
976:976.a:78262366554954104:[[1012,52318,-976,-37088]]:[[[1,2,0,-2,-1,-1],[0,0,1,1]]]:[29]:[[1012,52318,-976,-37088]]:[[[1,-2,0,2,-1],[0,0,1,1]]]:[[0]]
```

# Heuristically verifying the computation
It depends on [CHIMP](https://github.com/edgarcosta/CHIMP) and [smalljac](https://math.mit.edu/~drew/)
(no assert failing means success)
```
magma -b input:="704:704.a:2159575534599616393:[[-1856,129536,-44,-2948]]:[[[-1,-1,-2,1,1,-2],[1,1,1,1]]]:[2,3]:[[-1856,129536,-44,-2948],[771904,687149056,260876,-197548868]]:[[[0,4,4,-1,-2],[0,0,0,1]],[[-1140,-5166,-4646,2135,-236,4],[1]]]:[[0,9],[9,0]]" verify.m
```

