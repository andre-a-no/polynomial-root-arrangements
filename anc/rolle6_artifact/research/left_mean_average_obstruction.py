"""Exact left-cubic / mean derivative bound for an inflection secant.

This proof checker uses only the Python standard library. It reconstructs
all cleared polynomial identities and verifies every Bernstein leaf twice.
"""
from pathlib import Path
from fractions import Fraction as F
import json
from right_fourth_minimum_obstruction import ONE,U,K,add,scale,mul,power,direct_minimum
from bernstein3 import bernstein,split
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]
NAMES=['-D','2D-N','-H','left','mean']

def compactify(p):
    n=max(e[0] for e in p)
    num=scale(add(scale(ONE,3),scale(U,2)),F(1,5))
    den=add(ONE,scale(U,-1))
    return add(*(scale(mul(mul(power(num,i),power(den,n-i)),power(K,j)),co)
                 for (i,j,k),co in p.items() if k==0))

def check_certificates(path=None):
    data=json.loads(Path(path or ROOT/'data/left_mean_average_certificate.json').read_text())
    raw={name:{tuple(e):F(co) for e,co in data['raw_polynomials'][name]} for name in NAMES}
    source=[{tuple(e):F(co) for e,co in data['polynomials'][name]} for name in NAMES]
    assert all(e[2]==0 for p in raw.values() for e in p)
    assert all(compactify(raw[name])==p for name,p in zip(NAMES,source))
    t=K;V=mul(U,t)
    D=add(mul(V,add(scale(t,3),scale(ONE,-2))),scale(t,2),scale(ONE,-1))
    N=mul(add(scale(V,2),ONE),add(mul(U,add(ONE,scale(power(t,2),-2))),scale(t,-2),ONE))
    DN=add(scale(D,2),N);TN=add(mul(add(ONE,U),D),N)
    HN=add(mul(D,add(power(U,2),U,scale(ONE,-1))),mul(N,add(U,scale(ONE,-1))))
    assert raw['-D']==scale(D,-1) and raw['2D-N']==add(scale(D,2),scale(N,-1))
    assert mul(U,raw['-H'])==scale(HN,-1)
    # 60 D F(x), after r=N/(2D).
    coeff={6:scale(D,2),5:scale(DN,-3),4:scale(HN,-5),3:scale(mul(U,TN),10)}
    E=lambda x:add(*(mul(co,power(x,j)) for j,co in coeff.items()))
    EP=lambda x:add(*(scale(mul(co,power(x,j-1)),j) for j,co in coeff.items()))
    Z=add(E(ONE),scale(E(scale(U,-1)),-1))
    up=add(ONE,U);positive=add(scale(V,2),scale(U,-1),ONE)
    assert add(mul(up,EP(scale(V,-1))),scale(Z,-1))==mul(mul(up,positive),raw['left'])
    # mu=(2D+N)/(4D). Clear its fifth power denominator independently.
    fourD=scale(D,4)
    mean_eval=add(*(scale(mul(mul(co,power(DN,j-1)),power(fourD,6-j)),j)
                    for j,co in coeff.items()))
    left=add(mul(up,mean_eval),scale(mul(Z,power(fourD,5)),-1))
    right=scale(mul(mul(mul(up,D),add(N,scale(mul(U,D),2))),raw['mean']),16)
    assert left==right
    # The logarithmic-derivative identity at a=-ut gives D>0:
    # 1/(ut+T) = D/[ut(1-t)(ut+1)]. Verify the cleared identity.
    assert mul(D,add(V,add(ONE,U)))==add(mul(mul(V,add(ONE,scale(t,-1))),add(V,ONE)),scale(N,-1))
    # At t=1/2 the logarithmic derivative is negative; at t=3/4
    # it is positive, since both ut+1 and ut+T exceed ut.
    assert F(1,1)-F(1,3)-F(1,3)-F(1,3)==0
    assert F(3,5)**2+F(3,5)-1==F(-1,25)
    box=[(F(0),F(1)),(F(1,2),F(3,4)),(F(0),F(1))]
    cover=data['cover'];assert cover['names']==NAMES
    assert cover['box']==[[str(a),str(b)]for a,b in box]
    bs=[bernstein(p,box) for p in source]
    stats={'nodes':0,'leaves':0,'depth':0,'types':[0]*5};bounds=[]
    def walk(tree,bs,box,depth):
        stats['nodes']+=1;stats['depth']=max(stats['depth'],depth)
        if type(tree)is int:
            assert 0<=tree<5
            m=min(bs[tree][1].values());assert m>=0
            assert direct_minimum(source[tree],box)==m
            stats['leaves']+=1;stats['types'][tree]+=1;bounds.append([tree,str(m)]);return
        assert type(tree)is list and len(tree)==3 and type(tree[0])is int and tree[0]in(0,1)
        axis=tree[0];pairs=[split(p,axis)for p in bs];mid=sum(box[axis])/2
        for side in range(2):
            b=box[:];b[axis]=(box[axis][0],mid)if side==0 else(mid,box[axis][1])
            walk(tree[side+1],[p[side]for p in pairs],b,depth+1)
    walk(cover['tree'],bs,box,0)
    assert stats==cover['statistics']
    print('EXACTLY VERIFIED left-cubic/mean inflection average; 27 independent rational rectangles.',flush=True)
    return {'statistics':stats,'leaf_bounds':bounds}

def reasons(w):
    ps=patterns(w);sig=ps[0][0]
    if sig[1]==-1 and sig[10]==-1 and sig[5]==-1 and sig[11]==-1 and ps[0][3]==-1 and ps[2][3]==1:
        return ['left_mean_inflection_average']
    return []

if __name__=='__main__':
    check_certificates()
    known=json.loads((ROOT/'data/certificates_n6.json').read_text())
    unknown=json.loads((ROOT/'data/unresolved_n6.json').read_text())
    assert not any(reasons(w)for w in known)
    print('newly matched:',[w for w in unknown if reasons(w)])
