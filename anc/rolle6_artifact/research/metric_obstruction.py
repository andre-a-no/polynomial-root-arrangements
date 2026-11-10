"""Exact rational cubic-gap and derivative-span necessary inequalities.

The mathematical proof is in PROOF_METRIC_CONSTRAINTS.md. Certificate
verification uses Python integers and the standard library only.
"""
import json
from math import isqrt,gcd
from subfactor_obstruction import ROOT,FULL,SUBSETS,D,constraint_forms as base_forms
DIRECTIONS=[(1,1,56)]+[(1,k,64)for k in(2,4,8,16)]+[(k,1,64)for k in(2,4,8,16)]+[(3,5,64),(5,3,64)]
SPAN_CUTS={4:(70,99),5:(189,244),6:(218,267)}

def metric_rows(w):
    center=w.index('5');pos=[[i for i,c in enumerate(w)if c==str(j)]for j in range(6)]
    P=[[int(j<i)-int(j<center)for j in range(20)]+[0]*(D-20)for i in range(21)]
    Q={FULL:[[P[i]for i in row]for row in pos]};idx=20
    for S in SUBSETS:
        q=[[P[pos[0][i]]for i in S]]
        for j in range(1,len(S)):
            row=[]
            for i in range(len(S)-j):
                z=[0]*D;z[idx]=1;idx+=1;row.append(z)
            q.append(row)
        Q[S]=q
    assert idx==D
    out=[]
    def lin(terms):return [sum(c*v[i]for c,v in terms)for i in range(D)]
    for S in [FULL]+SUBSETS:
        if len(S)<3:continue
        j=len(S)-3;a,b,c=Q[S][j];u,v=Q[S][j+1]
        for t,s,m in DIRECTIONS:
            N=t*t+t*s+s*s;n=isqrt(m*m*N)+1;assert n*n>m*m*N
            den,left,right=3*n,m*(2*t+s),m*(t+2*s);g=gcd(den,gcd(left,right))
            den,left,right=den//g,left//g,right//g
            out.append(lin([(den,v),(-den,u),(-left,b),(left,a),(-right,c),(right,b)]))
    assert len(out)==462
    for S in [FULL]+SUBSETS:
        for j in range(max(0,len(S)-3)):
            m=len(S)-j;num,den=SPAN_CUTS[m];assert m*num*num<(m-2)*den*den
            x,y=Q[S][j],Q[S][j+1]
            out.append(lin([(den,y[-1]),(-den,y[0]),(-num,x[-1]),(num,x[0])]))
    assert len(out)==492
    return out

def constraint_forms(w,full=False):
    E,B=base_forms(w);rows=metric_rows(w)
    if not full:rows=rows[:462:11]+rows[462:]
    assert len(rows)==(492 if full else 72)
    return E,B+rows

def verify_certificate(w,record,full=False):
    E,B=constraint_forms(w,full=full);y=record['equalities'];z=record['positive']
    for pairs,size in [(y,len(E)),(z,len(B))]:
        assert len({i for i,a in pairs})==len(pairs)
        assert all(type(i)is int and type(a)is int and 0<=i<size and a!=0 for i,a in pairs)
    assert z and all(a>0 for i,a in z)
    for j in range(D):assert sum(a*E[i][j]for i,a in y)==sum(a*B[i][j]for i,a in z),(w,j)

def check_certificates():
    data=json.loads((ROOT/'data/metric_rejected_n6.json').read_text())
    for w,record in data.items():verify_certificate(w,record)
    alternate=json.loads((ROOT/'data/metric_full_rejected_n6.json').read_text())
    assert set(alternate)==set(data)
    for w,record in alternate.items():verify_certificate(w,record,full=True)
    rejected=set(data);assert {w[::-1]for w in rejected}==rejected
    print('EXACTLY VERIFIED',len(data),'metric obstruction certificates',flush=True)
    return rejected

if __name__=='__main__':check_certificates()
