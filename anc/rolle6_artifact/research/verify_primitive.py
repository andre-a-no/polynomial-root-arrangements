"""Independent exact checker for quintic-primitive obstruction identities.

Only Python's standard library is used. Polynomials in four gap variables
are dictionaries from exponent tuples to rational coefficients.
"""
import json,math
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ZERO=(0,0,0,0)

def constant(c):return {ZERO:Fraction(c)} if c else {}
def add(*polys):
    out={}
    for p in polys:
        for m,c in p.items():out[m]=out.get(m,0)+c
    return {m:c for m,c in out.items()if c}
def scale(p,k):return {m:c*k for m,c in p.items()if c*k}
def mul(p,q):
    out={}
    for u,a in p.items():
        for v,b in q.items():
            m=tuple(i+j for i,j in zip(u,v));out[m]=out.get(m,0)+a*b
    return {m:c for m,c in out.items()if c}
def power(p,k):
    out=constant(1)
    for _ in range(k):out=mul(out,p)
    return out
def primitive(p):
    den=math.lcm(*(c.denominator for c in p.values()))
    vals={m:int(c*den)for m,c in p.items()};g=math.gcd(*vals.values())
    return {m:Fraction(c//g)for m,c in vals.items()}
def exact_divide(p,q):
    rest=dict(p);out={};leading=max(q)
    while rest:
        m=max(rest);diff=tuple(a-b for a,b in zip(m,leading))
        assert min(diff)>=0,'non-polynomial division'
        term={diff:rest[m]/q[leading]};out=add(out,term)
        rest=add(rest,scale(mul(term,q),-1))
    return out
def evaluate(coeff,x):
    out={}
    for c in reversed(coeff):out=add(mul(out,x),c)
    return out
def derivative(coeff):return [scale(c,i)for i,c in enumerate(coeff)][1:]

X=[{tuple(int(i==j)for j in range(4)):Fraction(1)}for i in range(4)]
R=[{}]
for v in X:R.append(add(R[-1],v))
Q=[constant(1)]
for r in R:
    nextq=[{}for _ in range(len(Q)+1)]
    for i,p in enumerate(Q):nextq[i]=add(nextq[i],scale(mul(p,r),-1));nextq[i+1]=add(nextq[i+1],p)
    Q=nextq
F60=[{}]+[scale(p,Fraction(60,i+1))for i,p in enumerate(Q)]
H=[]
for lo,hi in ((0,3),(4,1)):
    diff=add(evaluate(F60,R[hi]),scale(evaluate(F60,R[lo]),-1))
    span=add(R[max(lo,hi)],scale(R[min(lo,hi)],-1))
    h=exact_divide(diff,power(span,3));assert mul(h,power(span,3))==diff
    H.append(h)
D=[Q]
for _ in range(4):D.append(derivative(D[-1]))
MEAN=scale(add(*R),Fraction(1,5))
shifted=[add(*(scale(mul(Q[j],power(MEAN,j-k)),math.comb(j,k))for j in range(k,6)))for k in range(6)]
AA,BB,CC=shifted[3],shifted[2],shifted[1]
RES=add(power(add(scale(CC,20),scale(power(AA,2),-5)),2),scale(mul(AA,power(BB,2)),160))
MEAN_HEIGHT=[]
for r in R:
    diff=add(evaluate(F60,MEAN),scale(evaluate(F60,r),-1))
    square=power(add(MEAN,scale(r,-1)),2)
    h=exact_divide(diff,square);assert mul(h,square)==diff
    MEAN_HEIGHT.append(h)

def constraints(w,extended=False,mean_sign=None):
    assert len(w)==15 and all(w.count(str(i))==5-i for i in range(5))
    pos=[[j for j,c in enumerate(w)if int(c)==i]for i in range(5)]
    polys=list(H)
    for k in range(2,5):
        for j,r in enumerate(R):
            s=(-1)**sum(z>pos[0][j]for z in pos[k])
            polys.append(primitive(scale(evaluate(D[k],r),s)))
    if extended:
        for k in (1,2):
            s=(-1)**sum(z>pos[4][0]for z in pos[k])
            polys.append(primitive(scale(evaluate(D[k],MEAN),s)))
        s=(-1)**sum(z>r for r in pos[1]for z in pos[3])
        polys.append(primitive(scale(RES,s)))
    if mean_sign is not None:
        assert mean_sign in(-1,1)and extended
        which=(0,2,4)if mean_sign==1 else(1,3)
        polys.extend(primitive(scale(MEAN_HEIGHT[i],mean_sign))for i in which)
    return polys

def check_file(name):
    data=json.loads((ROOT/'data'/name).read_text())
    for key,cert in data.items():
        w=key.split(':')[0]
        polys=constraints(w,cert.get('extended',False),cert.get('mean_sign'));out={};weight=0
        for coeff,ids,mon in cert['terms']:
            assert isinstance(coeff,int)and coeff>0
            assert len(mon)==4 and all(isinstance(i,int)and i>=0 for i in mon)
            if isinstance(ids,int):ids=[]if ids==-1 else[ids]
            assert all(isinstance(i,int)and 0<=i<len(polys)for i in ids)
            p=constant(1)
            for i in ids:p=mul(p,polys[i])
            if len(ids)>1:p=primitive(p)
            term=mul(p,{tuple(mon):Fraction(coeff)})
            assert all(sum(m)==cert['degree']for m in term)
            out=add(out,term);weight+=coeff
        assert weight>0 and not out,w
    reflected=set(data)|{(w.split(':')[0][::-1]+(':'+w.split(':')[1]if ':'in w else''))for w in data}
    print('EXACTLY VERIFIED',len(data),'primitive obstruction identities;',len(reflected),'oriented quintic words',flush=True)
    return reflected

def check():
    bad=check_file('primitive_cubic_rejected_n5.json')|check_file('primitive_extended_rejected_n5.json')
    five=set(json.loads((ROOT/'data/certificates_n5.json').read_text()))
    six=set(json.loads((ROOT/'data/certificates_n6.json').read_text()))
    tail=lambda w:''.join(str(int(c)-1)for c in w if c!='0')
    good=set(map(tail,six))
    assert bad<=five and not bad&good
    assert good|bad==five and len(good)==96 and len(bad)==20
    print('Quintic orders:',len(good),'have a certified sextic primitive;',len(bad),'cannot have one;',len(five-good-bad),'unresolved',flush=True)
    return bad

def check_mean():
    bad=check_file('primitive_mean_rejected_n5.json')
    six=set(json.loads((ROOT/'data/certificates_n6.json').read_text()))
    def signature(w):
        return ''.join(str(int(c)-1)for c in w if c!='0')+':'+str((-1)**w[:w.index('5')].count('0'))
    good=set(map(signature,six))
    tails={s.split(':')[0]for s in good}
    allpairs={w+':'+str(s)for w in tails for s in(-1,1)}
    assert len(good)==162 and len(bad)==30
    assert not good&bad and good|bad==allpairs
    print('EXACT CLASSIFICATION: 162 realizable (quintic order, primitive sign at mean) pairs; 30 additional pairs excluded',flush=True)
    return bad

if __name__=='__main__':check();check_mean()
