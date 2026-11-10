"""Standard-library-only constraints and exact sparse Farkas certificates."""
import itertools,json
from pathlib import Path
from verify_bounds import forms
from factor_deletion_obstruction import deletion_bounds
from subfactor_sign_rules import propagate
ROOT=Path(__file__).resolve().parents[1]
FULL=tuple(range(6))
SUBSETS=[s for k in range(5,1,-1)for s in itertools.combinations(FULL,k)]
D=20+sum(len(s)*(len(s)-1)//2 for s in SUBSETS)


def constraint_forms(w):
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
    e,b=forms(w,6);E=[v+[0]*(D-20)for v in e];B=[v+[0]*(D-20)for v in b]
    def linear(terms):return [sum(c*v[i]for c,v in terms)for i in range(D)]
    def diff(a,b):return [x-y for x,y in zip(a,b)]
    for S in SUBSETS:
        q=Q[S];n=len(S)
        for j in range(n-1):E.append(linear([(1,v)for v in q[j]]+[(-(n-j),q[-1][0])]))
        for row in q:B.extend(diff(y,x)for x,y in zip(row,row[1:]))
        for j in range(n-1):
            for i,z in enumerate(q[j+1]):B.extend([diff(z,q[j][i]),diff(q[j][i+1],z)])
        for j in range(n-2):
            m=n-j
            for k in range(1,m):
                B.append(linear([(m-k+1,q[j+1][k-1]),(-(m-k),q[j][k-1]),(-1,q[j][k])]))
                B.append(linear([(1,q[j][k-1]),(k,q[j][k]),(-(k+1),q[j+1][k-1])]))
            for k in range(1,m-1):
                B.append(linear([(k,v)for v in q[j][:k+1]]+[(-k-1,v)for v in q[j+1][:k]]))
                B.append(linear([(-k,v)for v in q[j][-(k+1):]]+[(k+1,v)for v in q[j+1][-k:]]))
        if n==5:
            removed=next(i for i in FULL if i not in S);bounds=deletion_bounds(w,removed);assert bounds is not None
            for j in range(1,5):
                for i,z in enumerate(q[j]):
                    lo,hi=bounds[j][i];B.extend([diff(z,P[lo]),diff(P[hi],z)])
    for S in [FULL]+SUBSETS:
        if len(S)<3:continue
        for deleted in S:
            T=tuple(i for i in S if i!=deleted)
            for j in range(1,len(T)):
                for i,z in enumerate(Q[T][j]):B.extend([diff(z,Q[S][j][i]),diff(Q[S][j][i+1],z)])
    for S in SUBSETS:
        for x in S:
            if x==5 or x+1 in S:continue
            T=tuple(sorted((set(S)-{x})|{x+1}))
            for j in range(1,len(S)):B.extend(diff(a,b)for a,b in zip(Q[T][j],Q[S][j]))
    info=propagate(w);assert 'conflict'not in info,(w,info)
    for S in SUBSETS:
        for j in range(1,len(S)):
            for i,(lo,hi)in enumerate(info['bounds'][S,j]):B.extend([diff(Q[S][j][i],P[lo]),diff(P[hi],Q[S][j][i])])
    return E,B


def verify_certificate(w,record):
    E,B=constraint_forms(w);y=record['equalities'];z=record['positive']
    for pairs,size in [(y,len(E)),(z,len(B))]:
        assert len({i for i,a in pairs})==len(pairs)
        assert all(type(i)is int and type(a)is int and 0<=i<size and a!=0 for i,a in pairs)
    assert z and all(a>0 for i,a in z)
    for j in range(D):assert sum(a*E[i][j]for i,a in y)==sum(a*B[i][j]for i,a in z),(w,j)


def check_certificates():
    data=json.loads((ROOT/'data/subfactor_rejected_n6.json').read_text())
    for w,record in data.items():verify_certificate(w,record)
    rejected=set(data);assert {w[::-1]for w in rejected}==rejected
    print('EXACTLY VERIFIED',len(data),'shared-subfactor sign obstruction certificates',flush=True)
    return rejected

if __name__=='__main__':check_certificates()
