"""Full marker-value comparisons at a stationary point, exact certificates."""
import json
from pathlib import Path
from fractions import Fraction as F
from stationary_obstruction import patterns,bounds,raw_tables as base_tables,BOX,dense,pat_key
from two_parameter_obstruction import add,scale,mul,power,dense_to_sparse,bernstein2,split2
from bernstein_combinations import elevate2
ROOT=Path(__file__).resolve().parents[1]

def raw_inequalities(pat):
    sig,iv,crit=pat;tb,outside=base_tables(iv)
    lo,hi=bounds(pat);one={(0,0):1};ineq={}
    for n in lo:ineq['L'+n]=(one,scale(tb[n],-1))
    for n in hi:ineq['H'+n]=(scale(one,-1),tb[n])
    K={(1,0):1};S={(0,1):1};D=add(power(K,2),scale(K,-1),one)
    AN=add(scale(one,2),scale(K,-2),scale(power(K,2),-1))
    CN=add(scale(one,-1),scale(K,-2),scale(power(K,2),2));DN=scale(add(AN,CN),-1)
    BN=add(scale(mul(CN,power(D,2)),3),scale(power(CN,3),-1))
    ys=[AN,scale(D,-1),CN,{},D,DN]
    if iv in(0,6):
        H=add(one,scale(S,-1));X=add(mul(H,ys[0 if iv==0 else 5]),scale(mul(S,D),-1 if iv==0 else 1))
        ys=[mul(H,Y)for Y in ys];D=mul(H,D);BN=mul(BN,power(H,3))
    else:X=add(ys[iv-1],mul(S,add(ys[iv],scale(ys[iv-1],-1))))
    Q60=add(scale(power(X,5),3),scale(mul(power(D,2),power(X,3)),-30),scale(mul(BN,power(X,2)),30))
    def P120(Y):return add(power(Y,6),scale(mul(power(D,2),power(Y,4)),-15),scale(mul(BN,power(Y,3)),20))
    for j,Y in enumerate(ys):
        for h,Z in enumerate(ys[j+1:],j+1):
            if sig[3*j]==sig[3*h]:continue
            sums=lambda n:add(*(mul(power(Y,n-i),power(Z,i))for i in range(n+1)))
            div=add(sums(5),scale(mul(power(D,2),sums(3)),-15),scale(mul(BN,sums(2)),20))
            assert mul(div,add(Y,scale(Z,-1)))==add(P120(Y),scale(P120(Z),-1))
            a=scale(add(Y,Z,scale(X,-2)),sig[3*h])
            b=scale(add(div,scale(Q60,-2)),sig[3*h])
            ineq['P'+str(j)+str(h)]=(a,b)
    return ineq,outside

def witness_polynomials(ineq,left,right):
    ai,bi=ineq[left];aj,bj=ineq[right]
    # Weights -aj, ai cancel the unknown coefficient Atilde.
    return ai,scale(aj,-1),add(mul(aj,bi),scale(mul(ai,bj),-1))

def weights_positive(a,b,box=BOX):
    if not all(v>=0 for c in(a,b)for row in c for v in row):return False
    N=max(len(a),len(b))-1;M=max(len(a[0]),len(b[0]))-1
    a=elevate2(a,N,M);b=elevate2(b,N,M)
    q=[[x+y for x,y in zip(r,s)]for r,s in zip(a,b)]
    # Nonzero weights are required also on internal subdivision boundaries.
    # The four boundary lines of the original rectangle are outside the
    # physical open parameter domain and need not have positive weights.
    if not any(v>0 for row in q for v in row):return False
    l,r,bot,top=box
    if l>BOX[0]and not any(v>0 for v in q[0]):return False
    if r<BOX[1]and not any(v>0 for v in q[-1]):return False
    if bot>BOX[2]and not any(row[0]>0 for row in q):return False
    if top<BOX[3]and not any(row[-1]>0 for row in q):return False
    for x,i in [(l,0),(r,-1)]:
        for y,j in [(bot,0),(top,-1)]:
            if BOX[0]<x<BOX[1]and BOX[2]<y<BOX[3]and q[i][j]<=0:return False
    return True

def witness_positive(bs,box=BOX):
    a,b,f=bs
    return all(v>=0 for row in f for v in row)and weights_positive(a,b,box)

def outside_combination(f,g):
    """Find exact lambda>=0 with f+lambda*g Bernstein-nonnegative."""
    N=max(len(f),len(g))-1;M=max(len(f[0]),len(g[0]))-1
    f=elevate2(f,N,M);g=elevate2(g,N,M);lo=F(0);hi=None
    for row,other in zip(f,g):
        for a,b in zip(row,other):
            if b>0:lo=max(lo,-a/b)
            elif b<0:hi=-a/b if hi is None else min(hi,-a/b)
            elif a<0:return None
    if hi is not None and lo>hi:return None
    assert lo>=0 and all(a+lo*b>=0 for row,other in zip(f,g)for a,b in zip(row,other))
    return lo

def check_combination_leaf(t,ws,box):
    assert set(t)=={'witness','outside','multiplier'}
    i=t['witness'];j=t['outside'];q=F(*t['multiplier'])
    assert type(i)is int and type(j)is int and 0<=i<len(ws)and 0<=j<len(ws)
    assert len(ws[i])==3 and len(ws[j])==1 and q>=0
    assert weights_positive(ws[i][0],ws[i][1],box)
    f,g=ws[i][2],ws[j][0];N=max(len(f),len(g))-1;M=max(len(f[0]),len(g[0]))-1
    f=elevate2(f,N,M);g=elevate2(g,N,M)
    assert all(a+q*b>=0 for row,other in zip(f,g)for a,b in zip(row,other))

def check_certificates():
    data=json.loads((ROOT/'data/full_stationary_certificates.json').read_text())
    assert data['rectangle']==[[1,1],[11,8],[0,1],[1,1]]
    total=0;seen=set()
    for cert in data['certificates']:
        pat=pat_key(cert);assert pat not in seen;seen.add(pat)
        ineq,outside=raw_inequalities(pat);initial=[]
        for r in cert['witnesses']:
            if r['pair']==['outside']:
                assert len(r['polynomials'])==1 and dense_to_sparse(r['polynomials'][0]['coefficients'])==outside
                assert not r['polynomials'][0]['factors']
            else:
                raw=witness_polynomials(ineq,*r['pair'])
                assert len(r['polynomials'])==3
                for idx,(p,pr)in enumerate(zip(raw,r['polynomials'])):
                    q=dense_to_sparse(pr['coefficients']);factor={(0,0):1}
                    if idx<2:assert not pr['factors']
                    for fr in pr['factors']:
                        fp=dense_to_sparse(fr['coefficients']);n=fr['exponent'];assert type(n)is int and n>0
                        if fr['kind']=='physical':assert fp==scale(outside,-1)
                        else:
                            assert fr['kind']=='bernstein';bs=bernstein2(fr['coefficients'],*BOX)
                            assert all(v>=0 for row in bs for v in row)and any(v>0 for row in bs for v in row)
                        factor=mul(factor,power(fp,n))
                    assert mul(factor,q)==p
            initial.append([bernstein2(p['coefficients'],*BOX)for p in r['polynomials']])
        stats={'nodes':0,'leaves':0,'maxdepth':0}
        def walk(t,ws,depth,box):
            stats['nodes']+=1;stats['maxdepth']=max(stats['maxdepth'],depth)
            if type(t)is int:
                assert 0<=t<len(ws)
                if len(ws[t])==1:assert all(v>=0 for row in ws[t][0]for v in row)
                else:assert witness_positive(ws[t],box)
                stats['leaves']+=1;return
            if type(t)is dict:
                check_combination_leaf(t,ws,box);stats['leaves']+=1;return
            assert type(t)is list and len(t)==3 and t[0]in(0,1)
            ps=[[split2(b,t[0])for b in w]for w in ws]
            l,r,bot,top=box;axis=t[0]
            left=(l,(l+r)/2,bot,top)if axis==0 else(l,r,bot,(bot+top)/2)
            right=((l+r)/2,r,bot,top)if axis==0 else(l,r,(bot+top)/2,top)
            walk(t[1],[[p[0]for p in w]for w in ps],depth+1,left)
            walk(t[2],[[p[1]for p in w]for w in ps],depth+1,right)
        walk(cert['tree'],initial,0,BOX);assert stats==cert['statistics'];total+=stats['leaves']
    print('EXACTLY VERIFIED',len(seen),'full stationary patterns and',total,'Bernstein rectangles',flush=True)

_patterns=None
def reasons(w):
    global _patterns
    if _patterns is None:
        data=json.loads((ROOT/'data/full_stationary_certificates.json').read_text())
        _patterns={pat_key(r):i for i,r in enumerate(data['certificates'])}
    return ['full_stationary_pattern_'+str(_patterns[p])for p in patterns(w)if p in _patterns]

if __name__=='__main__':
    check_certificates();known=json.loads((ROOT/'data/certificates_n6.json').read_text())
    assert not any(reasons(w)for w in known)
    print('EXACTLY VERIFIED: no known witness excluded',flush=True)
