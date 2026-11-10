"""Exact constraints and integer certificates for the mean-polar quartic.

All arithmetic in this module is standard-library integer arithmetic.
See PROOF_POLAR_MEAN.md for root intervals and positivity of the forms.
"""
import json
from pathlib import Path
from verify_bounds import forms
ROOT=Path(__file__).resolve().parents[1]


def root_bounds(w):
    pos=[[i for i,c in enumerate(w)if c==str(j)]for j in range(6)]
    a=pos[5][0];bounds=[]
    for j in range(4):
        r=pos[j];y=pos[j+1];row=[]
        for i in range(len(r)-1):
            if r[i+1]<a:row.append([y[i],r[i+1]])
            elif r[i]>a:row.append([r[i],y[i]])
        assert len(row)==4-j
        bounds.append(row)
    for j in range(4):
        n=4-j
        for x,c in enumerate(w):
            def sign(d):
                return 0 if int(c)==d else (-1)**sum(t>x for t in pos[d])
            s,t=sign(j),sign(j+1)
            u=t*((a>x)-(a<x))
            if s==0:known=-u
            elif u==0:known=-s
            elif s*u>0:known=-s
            else:continue
            assert known in (-1,1)
            left=sum(b[1]<=x for b in bounds[j])
            inside=[i for i,b in enumerate(bounds[j])if b[0]<x<b[1]]
            assert len(inside)<=1
            if inside:
                i=inside[0]
                if known==(-1)**(n-left):bounds[j][i][0]=x
                else:bounds[j][i][1]=x
            else:
                assert known==(-1)**(n-left),(w,j,x)
    return bounds


def constraint_forms(w):
    E0,B0=forms(w,6);center=w.index('5');d=30
    E=[v+[0]*10 for v in E0];B=[v+[0]*10 for v in B0]
    P=[[int(j<k)-int(j<center)for j in range(20)]+[0]*10 for k in range(21)]
    q=[];idx=20
    for j in range(4):
        q.append([])
        for i in range(4-j):
            q[-1].append([int(k==idx)for k in range(d)]);idx+=1
    def combine(terms):return [sum(c*v[k]for c,v in terms)for k in range(d)]
    for j in range(3):E.append(combine([(1,v)for v in q[j]]+[(-(4-j),q[3][0])]))
    bounds=root_bounds(w)
    for j in range(4):
        for i,z in enumerate(q[j]):
            lo,hi=bounds[j][i]
            B.append(combine([(1,z),(-1,P[lo])]))
            B.append(combine([(1,P[hi]),(-1,z)]))
        for x,y in zip(q[j],q[j][1:]):B.append(combine([(1,y),(-1,x)]))
    for j in range(3):
        for i,z in enumerate(q[j+1]):
            B.append(combine([(1,z),(-1,q[j][i])]))
            B.append(combine([(1,q[j][i+1]),(-1,z)]))
    for j in range(2):
        m=4-j
        for k in range(1,m):
            B.append(combine([(m-k+1,q[j+1][k-1]),(-(m-k),q[j][k-1]),(-1,q[j][k])]))
            B.append(combine([(1,q[j][k-1]),(k,q[j][k]),(-(k+1),q[j+1][k-1])]))
        for k in range(1,m-1):
            B.append(combine([(k,z)for z in q[j][:k+1]]+[(-k-1,z)for z in q[j+1][:k]]))
            B.append(combine([(-k,z)for z in q[j][-(k+1):]]+[(k+1,z)for z in q[j+1][-k:]]))
    return E,B


def verify_certificate(w,certificate):
    E,B=constraint_forms(w);v=certificate
    assert len(v)==len(E)+len(B) and all(type(x)is int for x in v)
    y=v[:len(E)];z=v[len(E):]
    assert min(z)>=0 and sum(z)>0
    for k in range(30):
        assert sum(a*row[k]for a,row in zip(y,E))==sum(a*row[k]for a,row in zip(z,B)),(w,k)


def check_certificates():
    data=json.loads((ROOT/'data/polar_mean_rejected_n6.json').read_text())
    for w,certificate in data.items():verify_certificate(w,certificate)
    rejected=set(data)
    assert {w[::-1]for w in rejected}==rejected
    print('EXACTLY VERIFIED',len(data),'mean-polar quartic obstruction certificates',flush=True)
    return rejected


if __name__=='__main__':check_certificates()
