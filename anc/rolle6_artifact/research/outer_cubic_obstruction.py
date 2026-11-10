"""Exact univariate obstruction using the outer and middle roots of P'''.
Only Python's standard library is used. See OUTER_CUBIC_PROOF.md.
"""
import json
from fractions import Fraction as F
from math import factorial
from pathlib import Path
from parametric_middle_obstruction import add,scale,mul,power,trim,bernstein
ROOT=Path(__file__).resolve().parents[1]
DATA=json.loads((ROOT/'data/outer_cubic_certificates.json').read_text())
D=[1,-1,1];AN=[2,-2,-1];CN=[-1,-2,2];DN=scale(add(AN,CN),-1)
BN=add(scale(mul(CN,power(D,2)),3),scale(power(CN,3),-1))


def forms():
    # (degree in z, parameter coordinate A/B/C/constant, numerator, D exponent)
    terms=[(6,3,[F(1,120)],0),(4,3,[F(-1,8)],0),
           (3,3,scale(BN,F(1,6)),3),(2,0,[F(1,2)],0),
           (1,1,[1],0),(0,2,[1],0)]
    out=[]
    for n in (AN,scale(D,-1),CN,[],D,DN):
        for j in range(3):
            row=[[]for _ in range(4)]
            for degree,coord,num,denpow in terms:
                if degree<j:continue
                p=degree-j
                v=scale(mul(mul(num,power(n,p)),power(D,6-denpow-p)),
                        120*factorial(degree)//factorial(p))
                row[coord]=add(row[coord],v)
            assert all(v.denominator==1 if isinstance(v,F) else isinstance(v,int)
                       for ps in row for v in ps)
            out.append(row)
    return out


def quotient(p,q):
    p=list(map(F,trim(p)));q=list(map(F,trim(q)))
    out=[F(0)]*max(0,len(p)-len(q)+1)
    while p and len(p)>=len(q):
        n=len(p)-len(q);v=p[-1]/q[-1];out[n]=v
        p=add(p,[0]*n+scale(q,-v))
    assert not p
    return trim(out)


def signature(w):
    pos3=[i for i,c in enumerate(w) if c=='3']
    if pos3[1]>w.index('5'):w=w[::-1]
    pos=[[i for i,c in enumerate(w)if c==str(j)]for j in range(6)]
    markers=sorted(pos[3]+pos[4]+pos[5])
    assert ''.join(w[i]for i in markers)=='343543'
    return tuple((-1)**sum(i>t for i in pos[j])for t in markers for j in range(3))


def reasons(w):
    ss=signature(w)
    return ['outer_cubic_B'] if all(a==0 or a==b for a,b in zip(DATA['signature'],ss)) else []


def check_certificates():
    # Cubic root identities: a+c+d=0, a^2+ac+c^2=3, a*c*d=-b.
    assert not add(AN,CN,DN)
    assert not add(power(AN,2),mul(AN,CN),power(CN,2),scale(power(D,2),-3))
    assert not add(mul(mul(AN,CN),DN),BN)
    fs=forms();certs={t['id']:t for t in DATA['identities']}
    assert len(certs)==len(DATA['identities'])
    mask=[0]*18
    for cert in certs.values():
        total=[[]for _ in range(4)]
        assert len({t['row']for t in cert['terms']})==len(cert['terms'])>0
        for term in cert['terms']:
            i=term['row'];sg=term['sign'];w=term['weight']
            assert 0<=i<18 and sg in (-1,1) and w and all(isinstance(a,int)for a in w)
            assert mask[i]in(0,sg);mask[i]=sg
            for j in range(4):total[j]=add(total[j],scale(mul(w,fs[i][j]),sg))
        assert total[:3]==[[],[],[]]
        assert not add(total[3],cert['residual'])
    assert mask==DATA['signature']
    guard=DATA['guard_polynomial'];assert guard==[-1,-5,5]
    end=F(1);applications=0
    for piece in DATA['pieces']:
        left,right=[F(*v)for v in piece['interval']]
        assert left==end and left<right
        branches=piece['branches']
        assert [t['guard_sign']for t in branches]in([0],[-1,1])
        for branch in branches:
            cert=certs[branch['identity']];g=branch['guard_sign'];strict=False
            for term in cert['terms']:
                w=term['weight'];bs=bernstein(w,left,right)
                if min(bs)<0:
                    assert g in(-1,1)
                    q=quotient(w,scale(guard,g))
                    assert min(bernstein(q,left,right))>=0
                elif min(bs)>0:strict=True
            assert strict # One coefficient is positive throughout the closed piece.
            assert min(bernstein(cert['residual'],left,right))>=0
            applications+=1
        end=right
    assert end==F(11,8)
    # (1+sqrt(3))/2<11/8, checked without a floating-point approximation.
    assert F(3)<(2*end-1)**2
    print('EXACTLY VERIFIED',len(certs),'outer-cubic identities and',applications,
          'interval branches covering the full normalized parameter domain',flush=True)

if __name__=='__main__':
    check_certificates()
    known=json.loads((ROOT/'data/certificates_n6.json').read_text())
    assert not any(reasons(w)for w in known)
    remaining=json.loads((ROOT/'data/unresolved_n6.json').read_text())
    assert not any(reasons(w)for w in remaining)
    rejected=json.loads((ROOT/'data/outer_cubic_rejected_n6.json').read_text())
    assert rejected and all(reasons(w)==v for w,v in rejected.items())
    print('EXACTLY VERIFIED',len(rejected),'stored additional exclusions',sorted(rejected))
