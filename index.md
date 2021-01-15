## UP IAM Center of Cryptography - _SAGE_ functions

Functions written in _SAGE_ for certain cryptographic topics. At the moment, the functions cover the following topics:


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
```python
def dualF(f,n):
    w=BooleanFunction(f).walsh_hadamard_transform();
    t=[];
    for x in w:
        if x==2^((n+1)/2):
            t.append(0);
        if x==-2^((n+1)/2):
            t.append(1);
    return t;
```
```python
def dualFshift(f,n,i):
    Sf=walsh_support(f,n);
    A=shift_space(Sf,Sf[i]);
    fd=dualF(f,n);
    A_s=sorted(A);
    tt=[fd[A.index(A_s[i])] for i in [0..len(A_s)-1]];
    return tt;
```
```python
def walsh_support(f,n):
    w=BooleanFunction(f).walsh_hadamard_transform();
    l=len(w);
    V=VectorSpace(GF(2),n);
    ws=[];
    for v in V:
        i=ZZ(list(v),base=2);
        if(w[i] !=0):
            ws.append(vector(GF(2),v));
    return ws;
```
```python
def shift_space(W,u):
    new=[];
    for v in W:
        new.append(vector(GF(2),v)+vector(GF(2),u));
    return new;
```
```python
def intersectionWS(ws,l,k):
     return set(ws[l[0]]).intersection(*[ws[l[i]] for i in [1..k-1]]);
```
```python
def intSize(ws,n,k):
        c2=Combinations(range(n),k).list();
        t=[];
        for x in c2:
            d=len(intersectionWS(tuple(x),k));
            if d not in t:
                t.append(d);
        return t;
```
```python
def inverseWHT(w,n):
    t=[];
    V=sorted(VectorSpace(GF(2),n));
    for u in V:
        s=0;
        for x in V:
            i=ZZ(list(x), base=2);
            s=s+w[i]*(-1)^(vector(u,GF(2))*vector(x,GF(2)));
        t.append((1-(2^(-n))*s)/2);
    return t;
```
```python
def newVBoolF(sF,n):
    Fnew=[];
    Fduals=[dualF(coordinates[i],n) for i in [0..(n-1)]]
    for i in [0..(n-1)]:
        Fnew.append(BooleanFunction(inverseWHT(semibent(sF[i],Fduals[i],n),n)));
    FnewList=[list(Fnew[i].truth_table(format='int')) for i in [0..(n-1)]]
    tt=[ZZ(list(transpose(matrix(GF(2),FnewList))[i]),base=2) for i in [0..len(FnewList[1])-1]]
    newF=SBox(tt);
    return newF;
 ```   
```python
def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3 
```
```python
def CCZ(S1,S2):
    n=len(S1);
    V=list(VectorSpace(GF(2),n));
    
    def toBin(a,n):
        return a.bits()+[0 for i in [0..n-len(a.bits())-1]];
    
    G1=[[1]+list(V[i])+toBin(S1[i],n) for i in [0..2^n-1]];
    G2=[[1]+list(V[i])+toBin(S2[i],n) for i in [0..2^n-1]];
    
    G1=matrix(GF(2),G1).transpose();
    G2=matrix(GF(2),G2).transpose();
    
    C1=LinearCode(G1);
    C2=LinearCode(G2);
    
    return C1.is_permutation_equivalent(C2);
```
```python
def is_vectorial_bent(F):
    n=len(F);
    l=[F[i].digits(base=2,padto=n) for i in [0..2^n-1]];
    l=list(transpose(matrix(GF(2),l)).echelon_form());
    s=span(l,GF(2));
    s=sorted(s);
    G=[BooleanFunction(list(s[i])).is_bent() for i in [1..2^((n/2))-1]];
    if ((False in G) or (len(s)!=2^(n/2))):
        return False;
    else:
        return True;
```
```python
def bent_dual(f,n):
    t=[];
    w=f.walsh_hadamard_transform();
    for x in w:
        if x==2^(n/2):
            t.append(0);
        else:
            t.append(1);
    return BooleanFunction(t);
```
```python
def Property_PTau(g,n):
    g1=bent_dual(g,n);
    FF=sorted(GF(2^n));
    l=[sorted([x,y]) for x in [1..2^n-1] for y in [x+1..2^n-1] if g1.derivative(x).derivative(y).algebraic_degree()==-1];
    G=Graph();
    G.add_edges(l);
    return G.cliques_maximum();
```
```python
def Property_PTau_dual(g1,n):
    FF=sorted(GF(2^n));
    l=[sorted([x,y]) for x in [1..2^n-1] for y in [x+1..2^n-1] if g1.derivative(x).derivative(y).algebraic_degree()==-1];
    G=Graph();
    G.add_edges(l);
    return G.cliques_maximum();

def reduced_var(x,U):
    return tuple((U[i]*x).trace() for i in [0..len(U)-1]);

def reduced_poly(fun,r):
    return BooleanFunction([fun(i) for i in r]);
```
```python
def stand_dual(f,n):
    t=[];
    w=BooleanFunction(f).walsh_hadamard_transform();
    for x in w:
        if x==0:
            t.append(0);
        else:
            t.append(1);
    return BooleanFunction(t);
```
```python
def is_subfield(K,F):
    bool=True;
    t1=[x in F for x in K];
    if False in t1:
        return False;
    t2=[];
    for x in K:
        for y in K:
            if ((x!=0) and (y!=0)):
                t2.append((x+y in K) and (x*y in K) and (-x in K) and (1/x in K));
    if False in t2:
        return False;
    else:
        return True;
```
```python
def PlateauedFromBent(sf,fd,n):
    f1=[(-1)^(fd[i]) for i in [0..len(fd)-1]];
    s1=[ZZ(list(x),base=2) for x in sf];
    t=[];
    j=0;
    for i in [0..2^n-1]:
        if i in s1:
            t.append(2^((n+1)/2)*f1[j]);
            j=j+1;
        else:
            t.append(0);
    return BooleanFunction(inverseWHT(t,n));
```
```python
def is_semibent(f,n):
    if((f.is_plateaued()) and (2^((n+1)/2) in f.walsh_hadamard_transform())):
        return True;
    else:
        return False;
```
```python
def walsh_spectrum(F):
    n=len(F);
    S=[];
    for i in [1..2^n-1]:
        f=F.component_function(i);
        S=S+list(f.walsh_hadamard_transform());
    K=sorted(list(set(S)));
    return [[K[i],S.count(K[i])] for i in [0..len(K)-1]];
```
```python
def is_in_M_strong(tt,F,n):
    s=[];
    for q in [1..2^6-1]:
        def der(L,q,a,b):
            return ttab(L.component_function(q).derivative(a).derivative(b));
        eqF=[[a,b] for a in [0..2^n-1] for b in [a+1..2^n-1] if der(tt,q,a,b)==der(F,q,a,b)];
        G=Graph();
        G.add_edges(eqF);
        C=G.cliques_maximal();
        cl=list(sage.graphs.cliquer.all_cliques(G, 2^(n/2), 2^(n/2)));
        V=VectorSpace(GF(2),n);
        cl1=[[sorted(V)[i] for i in s] for s in cl]
        b1=[V.subspace(s) for s in cl1];
        for K in b1:
            s.append(len(span(K))==2^(n/2));
    if True in s:
        return True;
    else:
        return False;

def is_in_M_weak(tt,F,n,q):
    s=[];
    def der(L,q,a,b):
        return ttab(L.component_function(q).derivative(a).derivative(b));
    eqF=[[a,b] for a in [0..2^n-1] for b in [a+1..2^n-1] if der(tt,q,a,b)==der(F,q,a,b)];
    G=Graph();
    G.add_edges(eqF);
    cl=list(sage.graphs.cliquer.all_cliques(G, 2^(n/2), 2^(n/2)));
    V=VectorSpace(GF(2),n);
    cl1=[[sorted(V)[i] for i in s] for s in cl]
    b1=[V.subspace(s) for s in cl1];
    for K in b1:
        s.append(len(span(K))==2^(n/2));
        if True in s:
            return True;
    return False;

def is_in_M_bool(f,n):
    s=[];
    for a in [0..2^n-1]:
        for b in [a+1..2^n-1]:
            if set(ttab(f.derivative(a).derivative(b)))=={0}:
                s.append([a,b]);
    G=Graph();
    G.add_edges(s);
    if G.clique_number()==2^n:
        return True;
    cl=list(sage.graphs.cliquer.all_cliques(G, 2^(n/2), 2^(n/2)));
    V=VectorSpace(GF(2),n);
    cl1=[[sorted(V)[i] for i in s] for s in cl]
    b1=[V.subspace(s) for s in cl1];
    for K in b1:
        s.append(len(span(K))==2^(n/2));
        if True in s:
            return True;
    return False;
```
```python
def cliqueSpace(f,n):
    s=[];
    for a in [0..2^n-1]:
        for b in [a+1..2^n-1]:
            if set(ttab(f.derivative(a).derivative(b)))=={0}:
                s.append([a,b]);
    G=Graph();
    G.add_edges(s);
    cl=list(sage.graphs.cliquer.all_cliques(G, 2^(n/2), 2^(n/2)));
    return cl;
```
```python
def GammaRank(F,n):
    Fn=sorted(GF(2^n));
    A=[];
    def id(F,Fn,a,b):
        return [[x+Fn[a],x+F[b]] for x in Fn];
    for x in [0..2^n-1]:
        for y in [0..2^n-1]:
            for a in [0..2^n-1]:
                for b in [0..2^n-1]:
                    if [Fn[x],Fn[y]] in id(F,Fn,a,b):
                        A.append(1);
                    else:
                        A.append(0);
    B=matrix(GF(2),2^(2*n),A);
    return B.rank();
```
```python
def ttab(f):
    return list(f.truth_table(format='int'))
```
```python
def isAB(M):
    T=list(transpose(matrix(M)));
    T1=[ZZ(list(x),base=2) for x in T];
    T1=SBox(T1);
    return T1.is_almost_bent();
```
```python
def phi(u,fd,gd,Sf,Sg):
    L=[s for s in Sf if s not in intersection(Sf,Sg)];
    tt=[];
    for x in L:
        i=Sf.index(x);
        a=fd[i];
        v=vector(GF(2),u)+vector(GF(2),x);
        j=Sg.index(tuple(v));
        b=gd[j];
        tt.append(b);
    return tt;

def M1(f2d,S2,L,S12):
    M11=[];
    M12=[];
    S2_i=[ZZ(list(x),base=2) for x in S2];
    S12_i=[ZZ(list(x),base=2) for x in S12];
    for u in sorted(VectorSpace(GF(2),5)):
        if u in S2:
            uL=shift_space(L,u);
            uL_i=[ZZ(list(x),base=2) for x in uL];
            if sorted(uL_i)==S12_i:
                k=[ZZ(list(x),base=2) for x in uL];
                if len(intersection(k,S2_i))>0:
                    M11.append(BooleanFunction([f2d[S2_i.index(i)] for i in k]).algebraic_normal_form())
        else:
            uO=shift_space(S12,u);
            uO_i=[ZZ(list(x),base=2) for x in uO];
            if sorted(uO_i)==S12_i:
                k=[ZZ(list(x),base=2) for x in uO];
                if len(intersection(k,S2_i))>0:
                    M12.append(BooleanFunction([f2d[S2_i.index(i)] for i in k]).algebraic_normal_form())
    return([M11,M12])

def M2(f2d,S2,L,S12,S2_bez_Om):
    M11=[];
    M12=[];
    S2_i=[ZZ(list(x),base=2) for x in S2];
    S12_i=[ZZ(list(x),base=2) for x in S12];
    S2_bez_Om_i=[ZZ(list(x),base=2) for x in S2_bez_Om];
    for u in sorted(VectorSpace(GF(2),5)):
        if u in S2:
            uL=shift_space(L,u);
            uL_i=[ZZ(list(x),base=2) for x in uL];
            if sorted(uL_i)==S2_bez_Om_i:
                print(uL_i);
                print();
                g=[f2d[S2_i.index(i)] for i in uL_i];
                print(g);
                M11.append(BooleanFunction(g).algebraic_normal_form())
        else:
            uO=shift_space(S12,u);
            uO_i=[ZZ(list(x),base=2) for x in uO];
            if sorted(uO_i)==S2_bez_Om_i:
                print(uL_i);
                print();
                g=[f2d[S2_i.index(i)] for i in uO_i];
                print(g);
                M12.append(BooleanFunction(g).algebraic_normal_form())
    return([M11,M12])
```
```python
def supp_eqn(Sf,n):
    X=var(['x'+str(i) for i in [0..n-1]]);
    eqn=[sum([X[i]*int(list(v)[i]) for i in [0..n-1]])==1 for v in Sf];
    return solve_mod(eqn,2);
```