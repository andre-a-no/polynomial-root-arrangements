"""Exact degree elevation and nonnegative combinations. Standard library only."""
from fractions import Fraction as F

def elevate1(v,N):
    v=list(v)
    while len(v)<=N:
        n=len(v)-1
        v=[v[0]]+[F(i,n+1)*v[i-1]+F(n+1-i,n+1)*v[i]for i in range(1,n+1)]+[v[-1]]
    return v

def elevate2(b,N,M):
    b=[elevate1(row,M)for row in b]
    return [list(row)for row in zip(*(elevate1(col,N)for col in zip(*b)))]

def positive_combination(bs,terms):
    assert terms and len({i for i,q in terms})==len(terms)
    assert all(type(i)is int and 0<=i<len(bs)and F(*q)>0 for i,q in terms)
    N=max(len(bs[i])-1 for i,q in terms);M=max(len(bs[i][0])-1 for i,q in terms)
    out=[[F(0)for j in range(M+1)]for i in range(N+1)]
    for i,q in terms:
        e=elevate2(bs[i],N,M);w=F(*q)
        out=[[v+w*z for v,z in zip(row,other)]for row,other in zip(out,e)]
    return all(v>=0 for row in out for v in row)

def used_indices(t):
    if type(t)is int:return {t}
    if type(t)is dict:
        assert set(t)=={'combination'}
        return {i for i,q in t['combination']}
    return used_indices(t[1])|used_indices(t[2])
