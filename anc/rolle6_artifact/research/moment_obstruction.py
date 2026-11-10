"""Exact certificates from second moments and their one-root deletion identity."""
import json
from math import gcd,isqrt
from functools import lru_cache
from metric_obstruction import ROOT,constraint_forms as base_forms
from root_subset_coordinates import FULL,SUBSETS,D,coordinates,linear,difference,verify_sparse
DIRECTIONS_DELETE=[(1,0)]+[(1,s)for s in(-8,-4,-2,-1,1,2,4,8)]+[(2,s)for s in(-1,1)]

def primitive(row):
    g=gcd(*row);assert g
    return [v//g for v in row]

def compositions(total,k):
    if k==1:yield(total,);return
    for a in range(total+1):
        for b in compositions(total-a,k-1):yield(a,)+b

@lru_cache(None)
def directions(m):
    out=[]
    for gaps in compositions(4,m-1):
        x=[0]
        for gap in gaps:x.append(x[-1]+gap)
        a=[m*v-sum(x)for v in x];g=gcd(*a);out.append(tuple(v//g for v in a))
    return out

def constraint_forms(w):
    E,B=base_forms(w);Q=coordinates(w);extra=[]
    for S in [FULL]+SUBSETS:
        for j in range(max(0,len(S)-2)):
            R=Q[S][j];m=len(R);delta=difference(Q[S][-2][1],Q[S][-2][0]);weights=[]
            for k in range(1,m):
                num=isqrt(1024**2*4*k*(m-k)//(m*m*(m-1)))+1
                assert num*num*m*m*(m-1)>1024**2*4*k*(m-k);weights.append(num)
            row=linear([(a,difference(y,x))for a,x,y in zip(weights,R,R[1:])]+[(-1024,delta)])
            extra.append(primitive(row))
            for a in directions(m):
                assert sum(a)==0 and any(a)
                N=m*(m-1)*sum(t*t for t in a);n=isqrt(64**2*N)+1;assert n*n>64**2*N
                extra.append(primitive(linear([(n,delta)]+[(-128*t,r)for t,r in zip(a,R)])))
    assert len(extra)==927
    for S in [FULL]+SUBSETS:
        n=len(S)
        if n<3:continue
        ds=difference(Q[S][-2][1],Q[S][-2][0]);mu=Q[S][-1][0]
        for i,deleted in enumerate(S):
            T=tuple(k for k in S if k!=deleted);dt=difference(Q[T][-2][1],Q[T][-2][0]);R=difference(Q[S][0][i],mu)
            M=n*(n-1)**2;alpha=(n-2)*(n-1)**2;beta=4*n
            for t,s in DIRECTIONS_DELETE:
                K=alpha*t*t+beta*s*s;num=isqrt(64**2*M*K)+1;assert num*num>64**2*M*K
                extra.append(primitive(linear([(num,ds),(-64*alpha*t,dt),(-64*beta*s,R)])))
            if i in(0,n-1):
                sign=-1 if i==0 else 1;num=isqrt(1024**2*(n-2)//n)+1;assert num*num*n>1024**2*(n-2)
                extra.append(primitive(linear([((n-1)*num,dt),(2048*sign,R),(-(n-1)*1024,ds)])))
    assert len(extra)==2727
    return E,B+extra

def verify_certificate(w,record):verify_sparse(*constraint_forms(w),record)

def check_certificates():
    data=json.loads((ROOT/'data/moment_rejected_n6.json').read_text())
    for w,record in data.items():verify_certificate(w,record)
    result=set(data);assert {w[::-1]for w in result}==result
    print('EXACTLY VERIFIED',len(result),'coupled second-moment certificates',flush=True)
    return result

if __name__=='__main__':check_certificates()
