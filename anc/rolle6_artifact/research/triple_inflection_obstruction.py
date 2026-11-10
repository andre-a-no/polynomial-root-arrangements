"""Exact obstruction involving three roots of P'' and the common mean.

Only standard-library rational arithmetic is used for the certificate.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import json
from two_parameter_obstruction import add,scale,mul,power,bernstein2,split2
from stationary_obstruction import dense
from inflection_obstruction import patterns

ROOT=Path(__file__).resolve().parents[1]
ONE={(0,0):F(1)}
U={(1,0):F(1)}
R={(0,1):F(1)}

def polynomials():
    H=add(power(U,2),U,scale(mul(U,R),2),scale(R,-2),scale(ONE,-1))
    G1={(0,4):5,(1,3):10,(2,2):-10,(1,2):-20,(0,2):-20,
        (3,1):20,(2,1):40,(1,1):-50,(0,1):-20,
        (4,0):24,(3,0):60,(2,0):10,(1,0):-20,(0,0):-5}
    G2={(0,4):16,(1,3):40,(0,3):16,(2,2):32,(1,2):24,
        (3,1):12,(2,1):14,(1,1):-6,(0,1):-4,
        (4,0):3,(3,0):6,(2,0):-1,(1,0):-4,(0,0):-1}
    return {'minus_H':scale(H,-1),'G1':G1,'G2':G2}

def identity_check():
    t=add(ONE,U,scale(R,2));m=scale(add(ONE,R),F(1,2))
    # Coefficients in x, with coefficients in Q[u,r].
    coeff=[ONE]
    for root in [scale(U,-1),{},ONE,t]:
        nxt=[{} for _ in range(len(coeff)+1)]
        for j,c in enumerate(coeff):
            nxt[j]=add(nxt[j],scale(mul(c,root),-1))
            nxt[j+1]=add(nxt[j+1],c)
        coeff=nxt
    primitive=[{},{}]+[scale(c,F(1,(j+1)*(j+2)))for j,c in enumerate(coeff)]
    def ev(x):return add(*(mul(c,power(x,j))for j,c in enumerate(primitive)))
    pu,pm,pt=ev(scale(U,-1)),ev(m),ev(t)
    ps=polynomials()
    lhs=scale(add(mul(m,pu),mul(U,pm)),960)
    rhs=scale(mul(mul(mul(U,m),add(R,scale(U,2),ONE)),ps['G1']),-1)
    assert lhs==rhs
    lhs=scale(add(mul(t,pu),mul(U,pt)),60)
    rhs=scale(mul(mul(mul(U,t),add(scale(R,2),scale(U,2),ONE)),ps['G2']),-1)
    assert lhs==rhs
    centered=[add(x,scale(m,-1))for x in [scale(U,-1),{},ONE,t]]
    assert not add(*centered)
    s2=scale(add(*(power(x,2)for x in centered)),F(1,12))
    assert scale(add(s2,scale(power(m,2),-1)),6)==scale(ps['minus_H'],-1)
    e3=add(*(mul(mul(centered[i],centered[j]),centered[k])
            for i in range(4)for j in range(i+1,4)for k in range(j+1,4)))
    assert e3==mul(mul(R,add(R,U)),add(R,U,ONE))
    # For u>=1, every coefficient of G2(1+h,r) is positive.
    shifted={}
    for (i,j),c in ps['G2'].items():
        for h in range(i+1):shifted[(h,j)]=shifted.get((h,j),0)+c*comb(i,h)
    shifted={ij:c for ij,c in shifted.items()if c}
    assert all(c>0 for c in shifted.values()) and shifted[(0,0)]==3
    return shifted

def initial():
    return [bernstein2(dense(p),F(0),F(1),F(0),F(1))for p in polynomials().values()]

def build_tree():
    def go(bs,depth=0):
        for i,p in enumerate(bs):
            if min(c for row in p for c in row)>=0:return i
        assert depth<30, 'No certificate found; not an exclusion'
        axis=depth%2
        pairs=[split2(p,axis)for p in bs]
        return [axis,go([p[0]for p in pairs],depth+1),go([p[1]for p in pairs],depth+1)]
    return go(initial())

def verify_tree(tree):
    stats={'nodes':0,'leaves':0,'maxdepth':0,'by_polynomial':{k:0 for k in polynomials()}}
    names=list(polynomials())
    def go(node,bs,depth=0):
        stats['nodes']+=1;stats['maxdepth']=max(stats['maxdepth'],depth)
        if type(node)is int:
            assert 0<=node<len(bs)
            assert min(c for row in bs[node]for c in row)>=0
            stats['leaves']+=1;stats['by_polynomial'][names[node]]+=1
            return
        assert type(node)is list and len(node)==3 and type(node[0])is int and node[0]in(0,1)
        axis,left,right=node;pairs=[split2(p,axis)for p in bs]
        go(left,[p[0]for p in pairs],depth+1);go(right,[p[1]for p in pairs],depth+1)
    go(tree,initial());return stats

def check_certificates():
    shifted=identity_check()
    d=json.loads((ROOT/'data/triple_inflection_certificate.json').read_text())
    assert d['rectangle']==[[0,1],[1,1],[0,1],[1,1]]
    expected={k:[[list(ij),int(v)]for ij,v in sorted(p.items())]for k,p in polynomials().items()}
    assert d['polynomials']==expected
    assert d['positive_tail']==[[list(ij),int(v)]for ij,v in sorted(shifted.items())]
    stats=verify_tree(d['tree']);assert d['statistics']==stats
    print('EXACTLY VERIFIED triple-inflection value identities, root geometry, positive unbounded tail, and',stats['leaves'],'rational Bernstein rectangles.',flush=True)
    return d

def reasons(w):
    ps=patterns(w);sig=ps[0][0]
    if sig[5]==-1 and sig[9]==1 and sig[11]==1 and ps[0][3]==1 and ps[1][3]==-1 and ps[3][3]==1:
        return ['three_inflections_and_mean']
    return []

if __name__=='__main__':
    import sys
    if '--write'in sys.argv:
        shifted=identity_check();tree=build_tree()
        d={'rectangle':[[0,1],[1,1],[0,1],[1,1]],
           'polynomials':{k:[[list(ij),int(v)]for ij,v in sorted(p.items())]for k,p in polynomials().items()},
           'positive_tail':[[list(ij),int(v)]for ij,v in sorted(shifted.items())],
           'tree':tree,'statistics':verify_tree(tree)}
        (ROOT/'data/triple_inflection_certificate.json').write_text(json.dumps(d,indent=2)+'\n')
    check_certificates()
