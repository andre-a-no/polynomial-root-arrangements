"""Alternating values at two inflections and the extreme zeros of P''' ."""
from pathlib import Path
from fractions import Fraction as F
import json
from right_fourth_minimum_obstruction import ONE,U,K,Y,add,scale,mul,power,translate,direct_minimum
from bernstein3 import bernstein,split
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]
NAMES=['-D','-N','N-2D','-H','-V','W-uV']

def check_certificates(path=None):
    data=json.loads(Path(path or ROOT/'data/outer_cubic_chord_certificate.json').read_text())
    parse=lambda rows:{tuple(e):F(co)for e,co in rows}
    source={name:parse(data['polynomials'][name])for name in NAMES}
    tail={name:parse(row)for name,row in data['tail'].items()}
    # First reconstruct the remainder of F-xF(1) modulo R' in variables u,r.
    r=K;s=add(ONE,r);T=add(ONE,U,scale(r,2))
    H=add(power(U,2),U,scale(mul(U,r),2),scale(r,-2),scale(ONE,-1))
    coeff={6:scale(ONE,F(1,30)),5:scale(s,F(-1,10)),4:scale(H,F(-1,12)),3:scale(mul(U,T),F(1,6))}
    F1=add(*coeff.values());coeff[1]=scale(F1,-1)
    # x^3=3(1+r)x^2/2+Hx/2-uT/4 at a zero of R'.
    pp=[(ONE,{},{})]
    for i in range(6):
        a,b,c=pp[-1]
        pp.append((scale(mul(c,mul(U,T)),F(-1,4)),add(a,scale(mul(c,H),F(1,2))),add(b,scale(mul(c,s),F(3,2)))))
    C,B,A=[add(*(mul(co,pp[i][j])for i,co in coeff.items()))for j in range(3)]
    h=add(scale(s,F(3,2)),scale(Y,-1))
    k=add(power(Y,2),scale(mul(s,Y),F(-3,2)),scale(H,F(-1,2)))
    V=add(B,mul(A,h));W=add(C,scale(mul(A,k),-1))
    Vc=scale(A,-1)
    assert tail['c_coefficient']==scale(Vc,960)
    half={}
    for(i,j,k0),co in V.items():
        e=(i,j,0);half[e]=half.get(e,F(0))+co*F(1,2)**k0
    half={e:co for e,co in half.items()if co}
    D0=add(scale(ONE,3),scale(r,6),scale(H,2))
    assert mul(r,tail['bound'])==add(scale(mul(D0,half),-1920),scale(mul(r,Vc),-960))
    for p in tail.values():
        q=translate(p,(F(2),F(0),F(0)))
        assert all(co>0 for co in q.values()) and q.get((0,0,0),0)>0
    # R'(1/2)=r/2, R'(mu)=-r(r+u)(r+u+1).
    cubic=lambda x:add(scale(power(x,3),4),scale(mul(s,power(x,2)),-6),scale(mul(H,x),-2),mul(U,T))
    assert cubic(scale(ONE,F(1,2)))==scale(r,F(1,2))
    assert cubic(scale(s,F(1,2)))==scale(mul(mul(r,add(r,U)),add(r,U,ONE)),-1)
    # Now replace r by N/(2D) and the middle cubic root by c=K.
    c=K;D=add(U,scale(mul(add(U,scale(ONE,-1)),c),-2),scale(power(c,2),-3))
    N=mul(add(scale(c,2),scale(ONE,-1)),add(power(U,2),U,scale(c,2),scale(power(c,2),-2)))
    HN=add(mul(D,add(power(U,2),U,scale(ONE,-1))),mul(N,add(U,scale(ONE,-1))))
    assert source['-D']==scale(D,-1)and source['-N']==scale(N,-1)
    assert source['N-2D']==add(N,scale(D,-2))and source['-H']==scale(HN,-1)
    def cleared(p):
        assert max(e[1]for e in p)<=5
        return scale(add(*(scale(mul(mul(mul(power(U,i),power(N,j)),power(scale(D,2),5-j)),power(K,k0)),co)
                           for(i,j,k0),co in p.items())),960)
    VV=cleared(V);WW=cleared(W)
    assert mul(add(scale(c,2),scale(ONE,-1)),source['-V'])==scale(VV,-1)
    assert source['W-uV']==add(WW,scale(mul(U,VV),-1))
    # D=0 at the physical substitution is impossible: N>0.
    assert add(scale(mul(D,power(c,3)),8),scale(mul(add(scale(D,2),N),power(c,2)),-6),scale(mul(HN,c),-4),scale(mul(U,add(mul(D,add(ONE,U)),N)),2))=={}
    box=[(F(3,5),F(2)),(F(1,2),F(3,5)),(F(0),F(1))]
    cover=data['cover'];assert cover['names']==NAMES and cover['box']==[[str(a),str(b)]for a,b in box]
    polys=list(source.values());bs=[bernstein(p,box)for p in polys]
    stats={'nodes':0,'leaves':0,'depth':0,'types':[0]*6}
    def walk(tree,bs,box,depth):
        stats['nodes']+=1;stats['depth']=max(stats['depth'],depth)
        if type(tree)is int:
            assert 0<=tree<6
            m=min(bs[tree][1].values());assert m>=0 and direct_minimum(polys[tree],box)==m
            stats['leaves']+=1;stats['types'][tree]+=1;return
        assert type(tree)is list and len(tree)==3 and type(tree[0])is int and tree[0]in(0,1)
        axis=tree[0];pairs=[split(p,axis)for p in bs];mid=sum(box[axis])/2
        for side in range(2):
            b=box[:];b[axis]=(box[axis][0],mid)if side==0 else(mid,box[axis][1])
            walk(tree[side+1],[p[side]for p in pairs],b,depth+1)
    walk(cover['tree'],bs,box,0);assert stats==cover['statistics']
    print('EXACTLY VERIFIED outer-cubic chord: six rational rectangles, independent coefficient checks and a positive tail.',flush=True)
    return stats

def reasons(w):
    p=patterns(w);sig=p[0][0]
    if sig[5]==-1 and sig[11]==1 and sig[0]==-1 and sig[15]==1 and p[1][3]==1 and p[2][3]==-1:
        return ['outer_cubic_chord']
    return []

if __name__=='__main__':
    check_certificates()
    known=json.loads((ROOT/'data/certificates_n6.json').read_text());unknown=json.loads((ROOT/'data/unresolved_n6.json').read_text())
    assert not any(reasons(w)for w in known)
    print('newly matched:',[w for w in unknown if reasons(w)])
