"""Two chord comparisons at the middle zero of P-triple-prime."""
from pathlib import Path
from fractions import Fraction as F
import json
from right_fourth_minimum_obstruction import ONE,U,K,add,scale,mul,power,translate,direct_minimum
from bernstein3 import bernstein,split
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]

def check_certificates(path=None):
 data=json.loads(Path(path or ROOT/'data/central_cubic_chord_certificate.json').read_text())
 source=[{tuple(e):F(v)for e,v in data['polynomials'][name]}for name in ['2D-N','-H','-G']]
 C=K;D=add(U,scale(mul(add(U,scale(ONE,-1)),C),-2),scale(power(C,2),-3))
 N=mul(add(scale(C,2),scale(ONE,-1)),add(power(U,2),U,scale(C,2),scale(power(C,2),-2)))
 A=add(scale(D,2),N);TN=add(mul(add(ONE,U),D),N)
 HN=add(mul(D,add(power(U,2),U,scale(ONE,-1))),mul(N,add(U,scale(ONE,-1))))
 assert source[0]==add(scale(D,2),scale(N,-1)) and source[1]==scale(HN,-1)
 assert add(scale(mul(D,power(C,3)),8),scale(mul(A,power(C,2)),-6),scale(mul(HN,C),-4),scale(mul(U,TN),2))=={}
 coeff={6:scale(D,2),5:scale(A,-3),4:scale(HN,-5),3:scale(mul(U,TN),10)}
 lhs=add(*(mul(co,add(scale(power(U,n-1),(-1)**n),power(C,n-1)))for n,co in coeff.items()))
 assert lhs==scale(mul(add(U,C),source[2]),-1)
 K1={(4,0,0):F(3),(3,0,0):F(10),(2,0,0):F(5),(1,0,0):F(-10),(0,0,0):F(-5)}
 assert all(v>0 for v in translate(K1,(F(1),F(0),F(0))).values())
 assert F(3,5)**2+F(3,5)-1==F(-1,25)
 assert F(2)*(1-F(3,5))-6*F(1,2)<0
 assert F(3,5)-2*(F(3,5)-1)*F(3,5)-3*F(3,5)**2==0
 box=[(F(3,5),F(1)),(F(1,2),F(3,5)),(F(0),F(1))]
 cover=data['cover'];assert cover['box']==[[str(a),str(b)]for a,b in box]
 assert cover['names']==['2D-N','-H','-G']
 bs=[bernstein(p,box)for p in source];stats={'nodes':0,'leaves':0,'depth':0,'types':[0,0,0]}
 def walk(tree,bs,box,depth):
  stats['nodes']+=1;stats['depth']=max(stats['depth'],depth)
  if type(tree)is int:
   assert 0<=tree<3;m=min(bs[tree][1].values());assert m>=0
   assert direct_minimum(source[tree],box)==m
   stats['leaves']+=1;stats['types'][tree]+=1;return
  assert type(tree)is list and len(tree)==3 and type(tree[0])is int and tree[0]in(0,1)
  axis=tree[0];pairs=[split(p,axis)for p in bs];mid=sum(box[axis])/2
  for side in range(2):
   b=box[:];b[axis]=(box[axis][0],mid)if side==0 else(mid,box[axis][1])
   walk(tree[side+1],[p[side]for p in pairs],b,depth+1)
 walk(cover['tree'],bs,box,0);assert stats==cover['statistics']
 print('EXACTLY VERIFIED central-cubic chord identities and',stats['leaves'],'independently checked rational rectangles.',flush=True)
 return stats

def reasons(w):
 ps=patterns(w);sig=ps[0][0]
 if sig[5]==-1 and sig[11]==-1 and sig[6]==1 and ps[0][3]==1 and ps[1][3]==-1:
  return ['central_cubic_chord']
 return []

if __name__=='__main__':
 check_certificates()
 known=json.loads((ROOT/'data/certificates_n6.json').read_text());unknown=json.loads((ROOT/'data/unresolved_n6.json').read_text())
 assert not any(reasons(w)for w in known)
 print('newly matched:',[w for w in unknown if reasons(w)])
