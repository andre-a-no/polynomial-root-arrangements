"""Integer-only verifier for coupled deletions and their common quartics."""
import json,itertools
from pathlib import Path
from verify_bounds import forms
from factor_deletion_obstruction import constraint_forms as single_forms
ROOT=Path(__file__).resolve().parents[1]


def constraint_forms(w,indices,monotone=False):
    indices=tuple(indices);assert 2<=len(indices)<=6
    assert list(indices)==sorted(set(indices))and all(type(k)is int and 0<=k<6 for k in indices)
    pairs=list(itertools.combinations(indices,2));d=20+10*len(indices)+6*len(pairs);center=w.index('5')
    P=[[int(j<i)-int(j<center)for j in range(20)]+[0]*(d-20)for i in range(21)]
    e,b=forms(w,6);E=[v+[0]*(d-20)for v in e];B=[v+[0]*(d-20)for v in b]
    positions=[i for i,c in enumerate(w)if c=='0'];G={};H={}
    def linear(terms):return [sum(c*v[k]for c,v in terms)for k in range(d)]
    def difference(x,y):return [a-b for a,b in zip(x,y)]
    for ii,k in enumerate(indices):
        start=20+10*ii;cols=list(range(20))+list(range(start,start+10));ee,bb=single_forms(w,k)
        for z in ee[len(e):]:
            row=[0]*d
            for col,val in zip(cols,z):row[col]=val
            E.append(row)
        for z in bb[len(b):]:
            row=[0]*d
            for col,val in zip(cols,z):row[col]=val
            B.append(row)
        q=[[P[i]for j,i in enumerate(positions)if j!=k]];idx=start
        for j in range(1,5):
            q.append([])
            for i in range(5-j):
                row=[0]*d;row[idx]=1;q[-1].append(row);idx+=1
        G[k]=q
    def add_chain(q):
        n=len(q)
        for j in range(n-1):E.append(linear([(1,v)for v in q[j]]+[(-(n-j),q[-1][0])]))
        for row in q:B.extend(difference(y,x)for x,y in zip(row,row[1:]))
        for j in range(n-1):
            for i,z in enumerate(q[j+1]):B.extend([difference(z,q[j][i]),difference(q[j][i+1],z)])
        for j in range(n-2):
            m=n-j
            for k in range(1,m):
                B.append(linear([(m-k+1,q[j+1][k-1]),(-(m-k),q[j][k-1]),(-1,q[j][k])]))
                B.append(linear([(1,q[j][k-1]),(k,q[j][k]),(-(k+1),q[j+1][k-1])]))
            for k in range(1,m-1):
                B.append(linear([(k,z)for z in q[j][:k+1]]+[(-k-1,z)for z in q[j+1][:k]]))
                B.append(linear([(-k,z)for z in q[j][-(k+1):]]+[(k+1,z)for z in q[j+1][-k:]]))
    for pp,(a,b)in enumerate(pairs):
        q=[[P[x]for k,x in enumerate(positions)if k not in(a,b)]];idx=20+10*len(indices)+6*pp
        for j in range(1,4):
            q.append([])
            for i in range(4-j):
                row=[0]*d;row[idx]=1;q[-1].append(row);idx+=1
        H[a,b]=q;add_chain(q)
        for parent in(a,b):
            for j in range(1,4):
                for i,z in enumerate(q[j]):B.extend([difference(z,G[parent][j][i]),difference(G[parent][j][i+1],z)])
    if monotone:
        for a,b in itertools.combinations(indices,2):
            for j in range(1,5):B.extend(difference(x,y)for x,y in zip(G[a][j],G[b][j]))
        for pair1,pair2 in itertools.combinations(pairs,2):
            keep1=[k for k in range(6)if k not in pair1];keep2=[k for k in range(6)if k not in pair2]
            if all(x>=y for x,y in zip(keep1,keep2)):
                for j in range(1,4):B.extend(difference(x,y)for x,y in zip(H[pair1][j],H[pair2][j]))
            elif all(x<=y for x,y in zip(keep1,keep2)):
                for j in range(1,4):B.extend(difference(x,y)for x,y in zip(H[pair2][j],H[pair1][j]))
    return E,B


def verify_certificate(w,record):
    E,B=constraint_forms(w,record['indices']);v=record['certificate'];d=len(E[0])
    assert len(v)==len(E)+len(B)and all(type(x)is int for x in v)
    y=v[:len(E)];z=v[len(E):]
    assert min(z)>=0 and sum(z)>0
    for j in range(d):
        assert sum(a*row[j]for a,row in zip(y,E))==sum(a*row[j]for a,row in zip(z,B)),(w,j)


def conditional_forms(w):
    E,B=forms(w,6);pos=[[i for i,c in enumerate(w)if c==str(j)]for j in range(6)];extra=[]
    if pos[0][2]<pos[2][1]:extra.append([(pos[0][1],2),(pos[0][3],1),(pos[0][4],1),(pos[0][2],-4)])
    if pos[0][3]>pos[2][2]:extra.append([(pos[0][3],4),(pos[0][1],-1),(pos[0][2],-1),(pos[0][4],-2)])
    for terms in extra:B.append([sum(c for i,c in terms if i>j)for j in range(20)])
    return E,B


def check_certificates():
    data=json.loads((ROOT/'data/joint_deletion_rejected_n6.json').read_text())
    for w,record in data.items():verify_certificate(w,record)
    rejected=set(data);assert {w[::-1]for w in rejected}==rejected
    direct=json.loads((ROOT/'data/joint_conditional_rejected_n6.json').read_text())
    assert set(direct)==rejected
    for w,v in direct.items():
        E,B=conditional_forms(w)
        assert len(v)==len(E)+len(B)and all(type(x)is int for x in v)
        y=v[:len(E)];z=v[len(E):]
        assert min(z)>=0 and sum(z)>0
        for j in range(20):
            assert sum(a*row[j]for a,row in zip(y,E))==sum(a*row[j]for a,row in zip(z,B)),(w,j)
    print('EXACTLY VERIFIED',len(data),'joint-deletion obstruction certificates',flush=True)
    return rejected

if __name__=='__main__':check_certificates()
