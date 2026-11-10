"""Two exact obstructions from reducing F and F' modulo F''''.

Only standard-library rational arithmetic is used for verification.
"""
from pathlib import Path
from fractions import Fraction as F
import json
from right_fourth_minimum_obstruction import ONE,U,K,add,scale,mul,power,translate,direct_minimum
from bernstein3 import bernstein
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]

def check_certificates(path=None):
    data=json.loads(Path(path or ROOT/'data/fourth_remainder_certificate.json').read_text())
    parse=lambda p:{n:{tuple(e):F(c) for e,c in row}for n,row in p.items()}
    av=parse(data['average']);ch=parse(data['chord'])
    r=K;s=add(ONE,r);T=add(ONE,U,scale(r,2))
    H=add(power(U,2),U,scale(mul(U,r),2),scale(r,-2),scale(ONE,-1))
    coeff={6:scale(ONE,F(1,30)),5:scale(s,F(-1,10)),4:scale(H,F(-1,12)),3:scale(mul(U,T),F(1,6))}
    F1=add(*coeff.values())
    # x^2=s*x+H/6 at a root of R''=12x^2-12s*x-2H.
    powers=[(ONE,{})]
    for i in range(6):
        b,a=powers[-1];powers.append((scale(mul(a,H),F(1,6)),add(b,mul(a,s))))
    BF=add(*(mul(co,powers[i][0])for i,co in coeff.items()))
    AF=add(*(mul(co,powers[i][1])for i,co in coeff.items()))
    BQ=add(*(scale(mul(co,powers[i-1][0]),i)for i,co in coeff.items()))
    AQ=add(*(scale(mul(co,powers[i-1][1]),i)for i,co in coeff.items()))
    assert av['A']==scale(AQ,20)
    assert av['C']==scale(add(F1,scale(BQ,-1)),120)
    assert mul(r,av['S'])==scale(add(F1,scale(mul(AQ,s),F(-1,2)),scale(BQ,-1)),60)
    delta=add(AF,scale(F1,-1))
    assert mul(r,ch['S'])==scale(delta,-180)
    assert ch['-H']==scale(H,-1)
    # J=(1+r)*3240*BF/H-540*delta; verify with H cleared.
    assert mul(H,ch['J'])==add(scale(mul(s,BF),3240),scale(mul(H,delta),-540))
    for p,shift in [(av['S'],F(7,10)),(ch['S'],F(2))]:
        q=translate(p,(shift,F(0),F(0)))
        assert all(c>0 for c in q.values())and q.get((0,0,0),0)>0
    small=[(F(3,5),F(7,10)),(F(0),F(1,3)),(F(0),F(1))]
    large=[(F(3,5),F(2)),(F(0),F(1)),(F(0),F(1))]
    bounds={}
    for name,p,box,expected in [('A',av['A'],small,F(882373,270000)),('C',av['C'],small,F(7347,625)),('J',ch['J'],large,F(32843,625))]:
        m=min(bernstein(p,box)[1].values());assert m==expected and m>0 and direct_minimum(p,box)==m
        bounds[name]=str(m)
    assert F(19,100)/(2*(1-F(7,10)))==F(19,60)<F(1,3)
    print('EXACTLY VERIFIED two fourth-derivative remainder rules: two coefficient tails and three positive polynomial bounds on two rectangles.',flush=True)
    return bounds

def average_reasons(w):
    p=patterns(w);sig=p[0][0]
    if sig[5]==-1 and sig[4]==1 and sig[13]==1 and p[1][3]==1 and p[2][3]==-1:
        return ['fourth_derivative_average']
    return []

def chord_reasons(w):
    p=patterns(w);sig=p[0][0]
    if sig[5]==-1 and sig[11]==1 and sig[3]==-1 and sig[12]==1 and p[1][3]==1 and p[2][3]==-1:
        return ['fourth_derivative_chord']
    return []

if __name__=='__main__':
    check_certificates()
    known=json.loads((ROOT/'data/certificates_n6.json').read_text())
    unknown=json.loads((ROOT/'data/unresolved_n6.json').read_text())
    for f in [average_reasons,chord_reasons]:
        assert not any(f(w)for w in known)
        print(f.__name__,[w for w in unknown if f(w)])
