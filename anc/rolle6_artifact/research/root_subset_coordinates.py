"""Exact coordinate forms for the already established shared-subfactor model."""
from subfactor_obstruction import FULL,SUBSETS,D

def coordinates(w):
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
    return Q

def linear(terms):return [sum(c*v[i]for c,v in terms)for i in range(D)]

def difference(x,y):return [a-b for a,b in zip(x,y)]

def verify_sparse(E,B,record):
    y,z=record['equalities'],record['positive']
    for pairs,size in [(y,len(E)),(z,len(B))]:
        assert len({i for i,a in pairs})==len(pairs)
        assert all(type(i)is int and type(a)is int and 0<=i<size and a!=0 for i,a in pairs)
    assert z and all(a>0 for i,a in z)
    for j in range(D):assert sum(a*E[i][j]for i,a in y)==sum(a*B[i][j]for i,a in z),j
