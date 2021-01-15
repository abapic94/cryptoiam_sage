## UP IAM Center of Cryptography - _SAGE_ functions

Functions written in _SAGE_ for certain cryptographic topics. At the moment, the functions cover the following topics:

```
- constructing semi-bent Boolean Functions via Walsh support and a dual bent function
- computing the dual of a vectorial Boolean function F whose coordinates are semi-bent
- computing the Walsh support of a Boolean Function
- finding the intersection of Walsh supports
- computing the inverse Walsh-Hadamard transform
- constructing new Vectorial Boolean functions via duals and Walsh supports
- BCT for non-permutations
- Determining if K is a subfield of F
- Determining if a vectorial Boolean function is bent
- Determining if a vectorial bent Boolean function is weakly/strongly outside M^#
- Determining the equation (if it exists) which defines a Walsh support of a given Boolean Function
```

```python
def semibent(S,ttD,n):
    ws=[0 for i in [0..(2^n)-1]];
    j=0;
    for x in S:
        i=ZZ(list(x), base=2);
        ws[i]=2^((n+1)/2)*((-1)^(int(ttD[j])));
        j=j+1;
    return ws;
```
