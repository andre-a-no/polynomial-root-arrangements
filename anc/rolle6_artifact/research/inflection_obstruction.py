"""Exact obstructions anchored at a zero of P''. Standard library only."""
import json
from pathlib import Path
from fractions import Fraction as F
from stationary_obstruction import BOX,dense
from two_parameter_obstruction import add,scale,mul,power,dense_to_sparse,bernstein2,split2
from bernstein_combinations import positive_combination,used_indices

ROOT=Path(__file__).resolve().parents[1]

def patterns(w):
    pos=[[i for i,c in enumerate(w)if c==str(j)]for j in range(6)]
    if pos[3][1]>pos[5][0]:
        w=w[::-1];pos=[[i for i,c in enumerate(w)if c==str(j)]for j in range(6)]
    markers=sorted(pos[3]+pos[4]+pos[5]);assert ''.join(w[x]for x in markers)=='343543'
    sign=lambda j,x:(-1)**sum(i>x for i in pos[j])
    sig=tuple(sign(j,x)for x in markers for j in range(3))
    return [(sig,sum(y<x for y in markers),1 if r%2==0 else -1,sign(0,x))for r,x in enumerate(pos[2])]

def bounds(pat):
    sig,iv,crit,pvalue=pat
    assert len(sig)==18 and all(v in(-1,1)for v in sig)
    assert iv in range(7)and crit in(-1,1)and pvalue in(-1,1)
    signs0=list(sig[::3])+[pvalue];signs1=list(sig[1::3])+[crit]
    lo=[];hi=[]
    for j,v in enumerate(signs1):(lo if v>0 else hi).append('W'+str(j))
    rank=lambda j:iv if j==6 else j+(j>=iv)
    for j in range(7):
        for h in range(j+1,7):
            if signs0[j]==signs0[h]:continue
            lower=signs0[j]*(rank(j)-rank(h))>0
            (lo if lower else hi).append('S'+str(j)+str(h))
    direct=['J'+str(j)+('p'if sig[3*j+2]>0 else 'm')for j in range(6)]
    return lo,hi,direct

def raw_tables(iv):
    assert iv in range(7)
    one={(0,0):1};K={(1,0):1};S={(0,1):1}
    D=add(power(K,2),scale(K,-1),one)
    AN=add(scale(one,2),scale(K,-2),scale(power(K,2),-1))
    CN=add(scale(one,-1),scale(K,-2),scale(power(K,2),2));DN=scale(add(AN,CN),-1)
    BN=add(scale(mul(CN,power(D,2)),3),scale(power(CN,3),-1))
    ys=[AN,scale(D,-1),CN,{},D,DN]
    if iv in(0,6):
        H=add(one,scale(S,-1))
        X=add(mul(H,ys[0 if iv==0 else 5]),scale(mul(S,D),-1 if iv==0 else 1))
        ys=[mul(H,Y)for Y in ys];BN=mul(BN,power(H,3));D=mul(H,D)
    else:X=add(ys[iv-1],mul(S,add(ys[iv],scale(ys[iv-1],-1))))
    def V(Y):return add(scale(power(Y,4),-15),scale(mul(power(D,2),power(Y,2)),90),scale(mul(BN,Y),-60))
    def Q60(Y):return add(scale(power(Y,5),3),scale(mul(power(D,2),power(Y,3)),-30),scale(mul(BN,power(Y,2)),30))
    def P120(Y):return add(power(Y,6),scale(mul(power(D,2),power(Y,4)),-15),scale(mul(BN,power(Y,3)),20))
    A=V(X);ys.append(X);tb={}
    for j,Y in enumerate(ys):
        tb['W'+str(j)]=add(scale(Q60(Y),-2),scale(mul(A,Y),-2))
        if j<6:
            tb['J'+str(j)+'p']=add(V(Y),scale(A,-1))
            tb['J'+str(j)+'m']=scale(tb['J'+str(j)+'p'],-1)
    for j,Y in enumerate(ys):
        for h,Z in enumerate(ys[j+1:],j+1):
            sums=lambda n:add(*(mul(power(Y,n-i),power(Z,i))for i in range(n+1)))
            div=add(sums(5),scale(mul(power(D,2),sums(3)),-15),scale(mul(BN,sums(2)),20))
            assert mul(div,add(Y,scale(Z,-1)))==add(P120(Y),scale(P120(Z),-1))
            tb['S'+str(j)+str(h)]=add(scale(div,-1),scale(mul(A,add(Y,Z)),-1))
    return tb,CN

def pat_key(r):return(tuple(r['signature']),r['interval'],r['critical_sign'],r['p_sign'])

def required_mask(cert):
    mask={}
    for idx in used_indices(cert['tree']):
        name=cert['comparisons'][idx]
        if name=='outside':continue
        if name.startswith('J'):
            j=int(name[1]);mask[3*j+2]=1 if name[2]=='p'else -1
            continue
        for typ,n in [(1,name.split('-')[0]),(-1,name.split('-')[1])]:
            if n[0]=='W':
                j=int(n[1])
                if j==6:assert cert['critical_sign']==typ
                else:mask[3*j+1]=typ
            else:
                assert n[0]=='S';j,h=map(int,n[1:]);iv=cert['interval']
                rank=lambda r:iv if r==6 else r+(r>=iv)
                sj=typ*(1 if rank(j)>rank(h)else -1)
                for r,v in [(j,sj),(h,-sj)]:
                    if r==6:assert cert['p_sign']==v
                    else:mask[3*r]=v
    assert all(cert['signature'][j]==v for j,v in mask.items())
    return tuple(sorted(mask.items()))

def check_certificates(path=None):
    data=json.loads(Path(path or ROOT/'data/inflection_certificates.json').read_text())
    assert data['rectangle']==[[1,1],[11,8],[0,1],[1,1]]
    cache={};initial={};seen=set();leaf_count=0
    for r in data['polynomials']:
        iv=r['interval'];name=r['id'];key=(iv,name);assert key not in initial
        if iv not in cache:cache[iv]=raw_tables(iv)
        tb,outside=cache[iv];p=dense_to_sparse(r['coefficients'])
        assert all(type(v)is int for row in r['coefficients']for v in row)
        if name=='outside':assert p==outside and not r['factors']
        else:
            raw=tb[name]if name.startswith('J')else add(tb[name.split('-')[0]],scale(tb[name.split('-')[1]],-1))
            factor={(0,0):1}
            for fr in r['factors']:
                fp=dense_to_sparse(fr['coefficients']);n=fr['exponent'];assert type(n)is int and n>0
                if fr['kind']=='physical':assert fp==scale(outside,-1)
                else:
                    assert fr['kind']=='bernstein'
                    bs=bernstein2(fr['coefficients'],*BOX)
                    assert all(v>=0 for row in bs for v in row)and any(v>0 for row in bs for v in row)
                factor=mul(factor,power(fp,n))
            assert mul(factor,p)==raw,key
        initial[key]=bernstein2(r['coefficients'],*BOX)
    for c in data['certificates']:
        pat=pat_key(c);assert pat not in seen;seen.add(pat)
        required_mask(c)
        lo,hi,direct=bounds(pat);names=c['comparisons'];assert len(names)==len(set(names))
        for name in names:
            if name=='outside':continue
            if name.startswith('J'):assert name in direct
            else:l,h=name.split('-');assert l in lo and h in hi
        stats={'nodes':0,'leaves':0,'maxdepth':0}
        def walk(t,bs,depth):
            stats['nodes']+=1;stats['maxdepth']=max(stats['maxdepth'],depth)
            if type(t)is int:
                assert 0<=t<len(bs)and all(v>=0 for row in bs[t]for v in row)
                stats['leaves']+=1;return
            if type(t)is dict:
                assert set(t)=={'combination'}and positive_combination(bs,t['combination'])
                stats['leaves']+=1;return
            assert type(t)is list and len(t)==3 and type(t[0])is int and t[0]in(0,1)
            pairs=[split2(b,t[0])for b in bs]
            walk(t[1],[p[0]for p in pairs],depth+1);walk(t[2],[p[1]for p in pairs],depth+1)
        walk(c['tree'],[initial[(c['interval'],name)]for name in names],0)
        assert c['statistics']==stats;leaf_count+=stats['leaves']
    print('EXACTLY VERIFIED',len(seen),'inflection patterns;',len(data['polynomials']),
          'polynomial identities;',leaf_count,'Bernstein rectangles',flush=True)
    return seen

_patterns=None
def reasons(w):
    global _patterns
    if _patterns is None:
        data=json.loads((ROOT/'data/inflection_certificates.json').read_text())
        _patterns=[(r['interval'],r['critical_sign'],r['p_sign'],required_mask(r))for r in data['certificates']]
    return ['inflection_pattern_'+str(idx)for sig,iv,crit,pval in patterns(w)
            for idx,(ri,rc,rp,mask)in enumerate(_patterns)
            if(iv,crit,pval)==(ri,rc,rp)and all(sig[j]==v for j,v in mask)]

if __name__=='__main__':
    check_certificates()
    known=json.loads((ROOT/'data/certificates_n6.json').read_text())
    assert not any(reasons(w)for w in known)
    print('EXACTLY VERIFIED: no known witness excluded',flush=True)
