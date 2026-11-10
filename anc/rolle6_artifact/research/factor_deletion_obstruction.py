"""Integer-only certificates from deleting one root of a hyperbolic sextic."""
import json
from pathlib import Path
from verify_bounds import forms
ROOT=Path(__file__).resolve().parents[1]

def deletion_bounds(w,k):
    pos=[[i for i,c in enumerate(w)if c==str(j)]for j in range(6)]
    r=pos[0][k];bounds=[[[x,x]for i,x in enumerate(pos[0])if i!=k]]
    for j in range(1,5):
        bounds.append([[a,b]for a,b in zip(pos[j],pos[j][1:])])
        n=5-j
        for x,c in enumerate(w):
            signp=lambda d:0 if int(c)==d else (-1)**sum(t>x for t in pos[d])
            if x==r:known=signp(j+1)
            else:
                prev=bounds[j-1]
                if any(a==b==x for a,b in prev):s=0
                elif any(a<x<b for a,b in prev):s=None
                else:s=(-1)**sum(a>=x if a!=b else a>x for a,b in prev)
                t=signp(j)
                if s is None:continue
                if t==0:known=-s
                elif s==0 or s*t<0:known=t
                else:continue
                known*=((x>r)-(x<r))
            assert known in(-1,1)
            left=sum(b<=x for a,b in bounds[j]);inside=[i for i,(a,b)in enumerate(bounds[j])if a<x<b]
            assert len(inside)<=1
            if inside:
                i=inside[0]
                if known==(-1)**(n-left):bounds[j][i][0]=x
                else:bounds[j][i][1]=x
            elif known!=(-1)**(n-left):return None
    return bounds


def constraint_forms(w,k):
    E0,B0=forms(w,6);center=w.index('5');d=30
    E=[v+[0]*10 for v in E0];B=[v+[0]*10 for v in B0]
    P=[[int(j<i)-int(j<center)for j in range(20)]+[0]*10 for i in range(21)]
    pos=[i for i,c in enumerate(w)if c=='0']
    q=[[P[i]for j,i in enumerate(pos)if j!=k]];idx=20
    for j in range(1,5):
        q.append([])
        for i in range(5-j):
            q[-1].append([int(h==idx)for h in range(d)]);idx+=1
    def combine(terms):return [sum(c*v[h]for c,v in terms)for h in range(d)]
    bounds=deletion_bounds(w,k);assert bounds is not None
    for j in range(4):E.append(combine([(1,v)for v in q[j]]+[(-(5-j),q[4][0])]))
    for j in range(5):
        for i,z in enumerate(q[j]):
            lo,hi=bounds[j][i]
            if lo!=hi:
                B.append(combine([(1,z),(-1,P[lo])]))
                B.append(combine([(1,P[hi]),(-1,z)]))
        for x,y in zip(q[j],q[j][1:]):B.append(combine([(1,y),(-1,x)]))
    for j in range(4):
        for i,z in enumerate(q[j+1]):
            B.append(combine([(1,z),(-1,q[j][i])]))
            B.append(combine([(1,q[j][i+1]),(-1,z)]))
    for j in range(3):
        m=5-j
        for l in range(1,m):
            B.append(combine([(m-l+1,q[j+1][l-1]),(-(m-l),q[j][l-1]),(-1,q[j][l])]))
            B.append(combine([(1,q[j][l-1]),(l,q[j][l]),(-(l+1),q[j+1][l-1])]))
        for l in range(1,m-1):
            B.append(combine([(l,z)for z in q[j][:l+1]]+[(-l-1,z)for z in q[j+1][:l]]))
            B.append(combine([(-l,z)for z in q[j][-(l+1):]]+[(l+1,z)for z in q[j+1][-l:]]))
    return E,B


def case_forms(w,k,word):
    E,B=constraint_forms(w,k);d=30;center=w.index('5')
    P=[[int(j<i)-int(j<center)for j in range(20)]+[0]*10 for i in range(21)]
    pos=[i for i,c in enumerate(w)if c=='0']
    q=[[P[i]for j,i in enumerate(pos)if j!=k]];idx=20
    for j in range(1,5):
        q.append([])
        for i in range(5-j):
            q[-1].append([int(h==idx)for h in range(d)]);idx+=1
    seen=[0]*5;seq=[]
    for c in word:
        j=int(c);seq.append(q[j][seen[j]]);seen[j]+=1
    for a,b in zip(seq,seq[1:]):B.append([y-x for x,y in zip(a,b)])
    return E,B


def verify_identity(w,E,B,v):
    assert len(v)==len(E)+len(B)and all(type(a)is int for a in v)
    y=v[:len(E)];z=v[len(E):]
    assert min(z)>=0 and sum(z)>0
    for j in range(30):
        assert sum(a*row[j]for a,row in zip(y,E))==sum(a*row[j]for a,row in zip(z,B)),(w,j)


def verify_certificate(w,record):
    k=record['root_index'];assert type(k)is int and 0<=k<6
    if 'certificate'in record:
        E,B=constraint_forms(w,k);verify_identity(w,E,B,record['certificate'])
    else:
        five=set(json.loads((ROOT/'data/certificates_n5.json').read_text()))
        assert len(five)==116 and set(record['quintic_cases'])==five
        for word,v in record['quintic_cases'].items():
            E,B=case_forms(w,k,word);verify_identity(w,E,B,v)


def check_certificates():
    data=json.loads((ROOT/'data/factor_deletion_rejected_n6.json').read_text())
    for w,record in data.items():verify_certificate(w,record)
    rejected=set(data)
    assert {w[::-1]for w in rejected}==rejected
    print('EXACTLY VERIFIED',len(data),'root-deletion obstruction certificates',flush=True)
    return rejected


if __name__=='__main__':check_certificates()
