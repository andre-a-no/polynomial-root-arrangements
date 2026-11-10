"""Independently checked stationary-value obstructions. Standard library only.

Candidate discovery is numerical, but this module reconstructs every identity
and verifies exact rational Bernstein covers; the discovery code is not trusted.
"""
import json
from pathlib import Path
from fractions import Fraction as F
from two_parameter_obstruction import add,scale,mul,power,dense_to_sparse,bernstein2,split2
from bernstein_combinations import positive_combination,used_indices

ROOT=Path(__file__).resolve().parents[1]
BOX=(F(1),F(11,8),F(0),F(1))

def patterns(w):
    pos=[[i for i,c in enumerate(w) if c==str(j)] for j in range(6)]
    if pos[3][1]>pos[5][0]:
        w=w[::-1]
        pos=[[i for i,c in enumerate(w) if c==str(j)] for j in range(6)]
    markers=sorted(pos[3]+pos[4]+pos[5])
    assert ''.join(w[x] for x in markers)=='343543'
    signs=tuple((-1)**sum(i>x for i in pos[j]) for x in markers for j in range(3))
    return [(signs,sum(y<x for y in markers),1 if r%2 else -1) for r,x in enumerate(pos[1])]

def bounds(pat):
    signs,interval,crit=pat
    assert len(signs)==18 and all(x in (-1,1) for x in signs)
    assert interval in range(7) and crit in (-1,1)
    lo=[];hi=[]
    for j in range(6):
        (lo if signs[3*j+1]*(1 if j>=interval else -1)>0 else hi).append('T'+str(j))
        (lo if signs[3*j+2]>0 else hi).append('V'+str(j))
        if signs[3*j]!=crit:
            (hi if crit>0 else lo).append('U'+str(j))
    (hi if crit>0 else lo).append('Vx')
    return lo,hi

def dense(p):
    n=max((i for i,j in p),default=0);m=max((j for i,j in p),default=0)
    return [[p.get((i,j),0) for j in range(m+1)] for i in range(n+1)]

def raw_tables(interval):
    assert interval in range(7)
    one={(0,0):1};K={(1,0):1};S={(0,1):1}
    D=add(power(K,2),scale(K,-1),one)
    AN=add(scale(one,2),scale(K,-2),scale(power(K,2),-1))
    CN=add(scale(one,-1),scale(K,-2),scale(power(K,2),2))
    DN=scale(add(AN,CN),-1)
    BN=add(scale(mul(CN,power(D,2)),3),scale(power(CN,3),-1))
    ys=[AN,scale(D,-1),CN,{},D,DN]
    if interval in (0,6):
        H=add(one,scale(S,-1))
        X=add(mul(H,ys[0 if interval==0 else 5]),scale(mul(S,D),-1 if interval==0 else 1))
        ys=[mul(H,Y)for Y in ys]
        BN=mul(BN,power(H,3));D=mul(H,D)
    else:
        X=add(ys[interval-1],mul(S,add(ys[interval],scale(ys[interval-1],-1))))
    def sums(Y,n,weighted=False):
        return add(*(scale(mul(power(Y,n-i),power(X,i)),i+1 if weighted else 1) for i in range(n+1)))
    def T(Y):
        return add(scale(sums(Y,4),-3),scale(mul(power(D,2),sums(Y,2)),30),scale(mul(BN,add(Y,X)),-30))
    def U(Y):
        return add(scale(sums(Y,4,True),-1),scale(mul(power(D,2),sums(Y,2,True)),15),scale(mul(BN,add(Y,scale(X,2))),-20))
    def V(Y):
        return add(scale(power(Y,4),-15),scale(mul(power(D,2),power(Y,2)),90),scale(mul(BN,Y),-60))
    def P120(Y):
        return add(power(Y,6),scale(mul(power(D,2),power(Y,4)),-15),scale(mul(BN,power(Y,3)),20))
    def Q60(Y):
        return add(scale(power(Y,5),3),scale(mul(power(D,2),power(Y,3)),-30),scale(mul(BN,power(Y,2)),30))
    table={}
    for j,Y in enumerate(ys):
        delta=add(Y,scale(X,-1))
        assert mul(T(Y),delta)==add(Q60(X),scale(Q60(Y),-1))
        assert mul(U(Y),power(delta,2))==add(scale(P120(Y),-1),P120(X),scale(mul(Q60(X),delta),2))
        table['T'+str(j)]=T(Y);table['U'+str(j)]=U(Y);table['V'+str(j)]=V(Y)
    table['Vx']=V(X)
    return table,CN

def pat_key(r):return (tuple(r['signature']),r['interval'],r['critical_sign'])

def required_mask(cert):
    """Only signs actually used by the proof, not all 18 discovery signs."""
    mask={};iv=cert['interval'];crit=cert['critical_sign']
    for idx in used_indices(cert['tree']):
        name=cert['comparisons'][idx]
        if name=='outside':continue
        for typ,b in [(1,name.split('-')[0]),(-1,name.split('-')[1])]:
            if b=='Vx':assert crit==-typ;continue
            j=int(b[1]);col={'T':1,'V':2,'U':0}[b[0]]
            sign=typ*(1 if j>=iv else -1)if b[0]=='T'else typ
            if b[0]=='U':assert typ==-crit
            key=3*j+col
            assert key not in mask or mask[key]==sign
            assert cert['signature'][key]==sign
            mask[key]=sign
    return tuple(sorted(mask.items()))

def check_certificates(path=None):
    path=path or ROOT/'data/stationary_certificates.json'
    data=json.loads(Path(path).read_text())
    assert data['rectangle']==[[1,1],[11,8],[0,1],[1,1]]
    cache={};initial={};patterns_seen=set();leaf_count=0
    for r in data['polynomials']:
        interval=r['interval'];name=r['id'];key=(interval,name)
        assert key not in cache
        if interval not in cache:cache[interval]=raw_tables(interval)
        table,outside=cache[interval]
        p=dense_to_sparse(r['coefficients'])
        assert all(type(x)is int for row in r['coefficients'] for x in row)
        if name=='outside':assert p==outside and not r['factors']
        else:
            lo,hi=name.split('-');raw=add(table[lo],scale(table[hi],-1))
            factor={(0,0):1}
            for fr in r['factors']:
                fp=dense_to_sparse(fr['coefficients']);exponent=fr['exponent']
                assert type(exponent)is int and exponent>0
                if fr['kind']=='physical':assert fp==scale(outside,-1)
                else:
                    assert fr['kind']=='bernstein'
                    bs=bernstein2(fr['coefficients'],*BOX)
                    assert all(x>=0 for row in bs for x in row) and any(x>0 for row in bs for x in row)
                factor=mul(factor,power(fp,exponent))
            assert mul(factor,p)==raw,(interval,name)
        cache[key]=r
        initial[key]=bernstein2(r['coefficients'],*BOX)
    for cert in data['certificates']:
        pat=pat_key(cert);assert pat not in patterns_seen;patterns_seen.add(pat)
        required_mask(cert)
        lo,hi=bounds(pat);interval=pat[1]
        names=cert['comparisons'];assert len(names)==len(set(names))
        for name in names:
            if name!='outside':
                l,h=name.split('-');assert l in lo and h in hi
        nodes=leaves=depthmax=0
        def walk(tree,bs,depth):
            nonlocal nodes,leaves,depthmax
            nodes+=1;depthmax=max(depthmax,depth)
            if type(tree)is int:
                assert 0<=tree<len(names)
                assert all(x>=0 for row in bs[tree] for x in row)
                leaves+=1;return
            if type(tree)is dict:
                assert set(tree)=={'combination'}and positive_combination(bs,tree['combination'])
                leaves+=1;return
            assert type(tree)is list and len(tree)==3
            axis,left,right=tree;assert type(axis)is int and axis in (0,1)
            pairs=[split2(b,axis)for b in bs]
            walk(left,[p[0]for p in pairs],depth+1)
            walk(right,[p[1]for p in pairs],depth+1)
        walk(cert['tree'],[initial[(interval,n)] for n in names],0)
        assert cert['statistics']=={'nodes':nodes,'leaves':leaves,'maxdepth':depthmax}
        leaf_count+=leaves
    print('EXACTLY VERIFIED',len(patterns_seen),'stationary patterns;',len(data['polynomials']),
          'polynomial identities;',leaf_count,'Bernstein rectangles',flush=True)
    return patterns_seen

_patterns=None
def reasons(w):
    global _patterns
    if _patterns is None:
        data=json.loads((ROOT/'data/stationary_certificates.json').read_text())
        _patterns=[(r['interval'],r['critical_sign'],required_mask(r))for r in data['certificates']]
    return ['stationary_pattern_'+str(idx)for sig,iv,crit in patterns(w)
            for idx,(ri,rc,mask)in enumerate(_patterns)
            if iv==ri and crit==rc and all(sig[j]==v for j,v in mask)]

if __name__=='__main__':
    check_certificates()
    known=json.loads((ROOT/'data/certificates_n6.json').read_text())
    assert not any(reasons(w)for w in known)
    print('EXACTLY VERIFIED: no known witness excluded',flush=True)
