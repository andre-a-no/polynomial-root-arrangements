"""Exact obstruction at a positive stationary value between c and the mean.
Five rational comparisons, two parameter charts combined by a blow-up,
and an exact Bernstein subdivision certificate. Standard library only.
See PROOF_TWO_PARAMETERS.md.
"""
import json
from fractions import Fraction as F
from math import comb
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DATA=json.loads((ROOT/'data/two_parameter_certificates.json').read_text())

# Sparse bivariate polynomial arithmetic: (power of k, power of w) -> integer.
def add(*ps):
    r={}
    for p in ps:
        for ij,a in p.items():r[ij]=r.get(ij,0)+a
    return {ij:a for ij,a in r.items()if a}
def scale(p,a):return {ij:a*b for ij,b in p.items()if a*b}
def mul(p,q):
    r={}
    for (i,j),a in p.items():
        for (h,l),b in q.items():r[i+h,j+l]=r.get((i+h,j+l),0)+a*b
    return {ij:a for ij,a in r.items()if a}
def power(p,n):
    r={(0,0):1}
    for _ in range(n):r=mul(r,p)
    return r
def dense_to_sparse(a):return {(i,j):v for i,row in enumerate(a)for j,v in enumerate(row)if v}

def generated_comparisons():
    one={(0,0):1};K={(1,0):1};W={(0,1):1}
    km=add(K,scale(one,-1));D=add(power(K,2),scale(K,-1),one)
    AN=add(scale(one,2),scale(K,-2),scale(power(K,2),-1))
    CN=add(scale(one,-1),scale(K,-2),scale(power(K,2),2))
    BN=add(scale(mul(CN,power(D,2)),3),scale(power(CN,3),-1))
    G=mul(W,D);XN=add(scale(mul(K,km),3),scale(G,-1))
    UN=scale(G,-1);AN=mul(AN,W);CN=mul(CN,W);BN=mul(BN,power(W,3))
    # Every form below is multiplied by the positive denominator 60*(wD)^4.
    def sums(Y,n,weighted=False):
        return add(*(scale(mul(power(Y,n-i),power(XN,i)),i+1 if weighted else 1)for i in range(n+1)))
    def T(Y):
        return add(scale(sums(Y,4),-3),scale(mul(power(G,2),sums(Y,2)),30),scale(mul(BN,add(Y,XN)),-30))
    def U(Y):
        return add(scale(sums(Y,4,True),-1),scale(mul(power(G,2),sums(Y,2,True)),15),scale(mul(BN,add(Y,scale(XN,2))),-20))
    def V(Y):
        return add(scale(power(Y,4),-15),scale(mul(power(Y,2),power(G,2)),90),scale(mul(BN,Y),-60))
    fs={'Tu-Uc':add(T(UN),scale(U(CN),-1)),
        'Tc-Vu':add(T(CN),scale(V(UN),-1)),
        'Tc-Uu':add(T(CN),scale(U(UN),-1)),
        'Tu-Uu':add(T(UN),scale(U(UN),-1)),
        'Tc-Ua':add(T(CN),scale(U(AN),-1))}
    common=scale(power(km,3),3)
    factors={name:mul(common,power(K,3))for name in fs};factors['Tc-Ua']=common
    return XN,fs,factors


def bernstein2(p,left,right,bottom,top):
    n=len(p)-1;m=len(p[0])-1
    a=[[sum(p[I][J]*comb(I,i)*left**(I-i)*(right-left)**i*
          comb(J,j)*bottom**(J-j)*(top-bottom)**j
          for I in range(i,n+1)for J in range(j,m+1))for j in range(m+1)]for i in range(n+1)]
    return [[sum(a[i][j]*F(comb(I,i),comb(n,i))*F(comb(J,j),comb(m,j))
             for i in range(I+1)for j in range(J+1))for J in range(m+1)]for I in range(n+1)]

def split1(row):
    levels=[row]
    while len(levels[-1])>1:
        q=levels[-1];levels.append([(a+b)/2 for a,b in zip(q,q[1:])])
    return [v[0]for v in levels],[v[-1]for v in levels[::-1]]

def split2(p,axis):
    if axis==0:
        pairs=[split1(list(col))for col in zip(*p)]
        return [list(row)for row in zip(*(a for a,b in pairs))],[list(row)for row in zip(*(b for a,b in pairs))]
    pairs=[split1(row)for row in p]
    return [a for a,b in pairs],[b for a,b in pairs]

def reasons(word):
    w=word
    if [i for i,c in enumerate(w)if c=='3'][1]>w.index('5'):w=w[::-1]
    pos=[[i for i,c in enumerate(w)if c==str(j)]for j in range(6)]
    a,c,d=pos[3];u=pos[4][0];mu=pos[5][0]
    assert a<u<c<mu<pos[4][1]<d
    def sign(j,t):return (-1)**sum(i>t for i in pos[j])
    if any(sign(j,t)!=-1 for j,t in ((0,a),(0,u),(0,c),(1,u),(1,c),(2,u))):return []
    if not any(c<x<mu and sign(0,x)>0 for x in pos[1]):return []
    return ['positive_stationary_value']

def check_certificates():
    records=DATA['polynomials'];names=[r['id']for r in records]
    assert len(set(names))==len(names)
    polys={r['id']:dense_to_sparse(r['coefficients'])for r in records}
    assert all(isinstance(x,int)for r in records for row in r['coefficients']for x in row)
    outside,fs,factors=generated_comparisons()
    assert polys['outside']==outside
    assert set(polys)=={'outside','corner_combination'}|set(fs)
    for name,f in fs.items():assert mul(factors[name],polys[name])==f,name
    rec=next(r for r in records if r['id']=='corner_combination')
    assert rec['combination']=={'outside':45,'Tc-Vu':1}
    assert polys['corner_combination']==add(scale(outside,45),polys['Tc-Vu'])
    assert DATA['rectangle']==[[1,1],[11,8],[0,1],[1,1]]
    rect=[F(*v)for v in DATA['rectangle']]
    assert 3<(2*rect[1]-1)**2
    initial=[bernstein2(r['coefficients'],*rect)for r in records]
    counts={name:0 for name in names};nodes=0;maxdepth=0
    def walk(tree,bs,depth):
        nonlocal nodes,maxdepth
        nodes+=1;maxdepth=max(depth,maxdepth)
        if type(tree)is int:
            assert 0<=tree<len(records)
            assert min(x for row in bs[tree]for x in row)>=0
            counts[names[tree]]+=1
            return
        assert isinstance(tree,list)and len(tree)==3
        axis,left,right=tree;assert type(axis)is int and axis in (0,1)
        pairs=[split2(b,axis)for b in bs]
        walk(left,[v[0]for v in pairs],depth+1)
        walk(right,[v[1]for v in pairs],depth+1)
    walk(DATA['tree'],initial,0)
    stats={'nodes':nodes,'leaves':sum(counts.values()),'maxdepth':maxdepth,'by_polynomial':counts}
    assert stats==DATA['statistics']
    print('EXACTLY VERIFIED 5 two-parameter identities, one positive combination, and',
          stats['leaves'],'Bernstein rectangles covering the full domain',flush=True)

if __name__=='__main__':
    check_certificates()
    known=json.loads((ROOT/'data/certificates_n6.json').read_text())
    assert not any(reasons(w)for w in known)
    remaining=json.loads((ROOT/'data/unresolved_n6.json').read_text())
    assert not any(reasons(w)for w in remaining)
    rejected=json.loads((ROOT/'data/two_parameter_rejected_n6.json').read_text())
    assert rejected and all(reasons(w)==v for w,v in rejected.items())
    print('EXACTLY VERIFIED',len(rejected),'stored additional exclusions',sorted(rejected))
