"""Chord through three inflections and the first zero of P-triple-prime."""
from pathlib import Path
from fractions import Fraction as F
import json
from right_fourth_minimum_obstruction import ONE,U,K,add,scale,mul,power,translate,direct_minimum
from bernstein3 import bernstein,split
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]

def check_certificates(path=None):
 data=json.loads(Path(path or ROOT/'data/left_cubic_chord_certificate.json').read_text())
 source=[{tuple(e):F(v)for e,v in data['polynomials'][name]}for name in ['2D-N','-H','-G']]
 t=K;V=mul(U,t)
 D=add(mul(mul(U,t),add(scale(t,3),scale(ONE,-2))),scale(t,2),scale(ONE,-1))
 small=add(mul(U,add(ONE,scale(power(t,2),-2))),scale(t,-2),ONE)
 N=mul(add(scale(V,2),ONE),small)
 assert add(mul(add(ONE,scale(power(t,2),-2)),D),mul(mul(t,add(scale(ONE,2),scale(t,-3))),small))==mul(add(scale(t,2),scale(ONE,-1)),power(add(ONE,scale(t,-1)),2))
 A=add(scale(D,2),N);TN=add(mul(add(ONE,U),D),N)
 HN=add(mul(D,add(power(U,2),U,scale(ONE,-1))),mul(N,add(U,scale(ONE,-1))))
 assert source[0]==add(scale(D,2),scale(N,-1)) and source[1]==scale(HN,-1)
 # The first cubic root is x=-ut.
 assert add(scale(mul(D,power(scale(V,-1),3)),8),scale(mul(A,power(V,2)),-6),scale(mul(HN,V),4),scale(mul(U,TN),2))=={}
 coeff={6:scale(D,2),5:scale(A,-3),4:scale(HN,-5),3:scale(mul(U,TN),10)}
 E=lambda x:add(*(mul(co,power(x,n))for n,co in coeff.items()))
 chord=add(mul(add(ONE,V),E(scale(U,-1))),mul(add(U,scale(V,-1)),E(ONE)),scale(mul(add(ONE,U),E(scale(V,-1))),-1))
 assert chord==scale(source[2],-1)
 # Endpoint evaluations proving 1/2<t<2/3. r is represented by the third variable.
 R={(0,0,1):1}
 primitive_cubic=mul(U,add(N,scale(mul(R,D),-2)))
 def at_t(p,q):
  out={}
  for(i,j,l),v in p.items():out[(i,0,l)]=out.get((i,0,l),0)+v*q**j
  return {e:v for e,v in out.items()if v}
 assert at_t(primitive_cubic,F(1,2))==scale(mul(power(U,2),add(U,ONE,R)),F(1,2))
 assert at_t(primitive_cubic,F(2,3))==mul(U,add(scale(mul(add(scale(U,4),scale(ONE,3)),add(U,scale(ONE,-3))),F(1,27)),scale(R,F(-2,3))))
 box=[(F(3,5),F(1)),(F(1,2),F(2,3)),(F(0),F(1))]
 cover=data['cover'];assert cover['box']==[[str(a),str(b)]for a,b in box]
 assert cover['names']==['2D-N','-H','-G']
 bs=[bernstein(p,box)for p in source];stats={'nodes':0,'leaves':0,'depth':0,'types':[0,0,0]}
 def walk(tree,bs,box,depth):
  stats['nodes']+=1;stats['depth']=max(stats['depth'],depth)
  if type(tree)is int:
   assert 0<=tree<3;m=min(bs[tree][1].values());assert m>=0 and direct_minimum(source[tree],box)==m
   stats['leaves']+=1;stats['types'][tree]+=1;return
  assert type(tree)is list and len(tree)==3 and type(tree[0])is int and tree[0]in(0,1)
  axis=tree[0];pairs=[split(p,axis)for p in bs];mid=sum(box[axis])/2
  for side in range(2):
   b=box[:];b[axis]=(box[axis][0],mid)if side==0 else(mid,box[axis][1])
   walk(tree[side+1],[p[side]for p in pairs],b,depth+1)
 walk(cover['tree'],bs,box,0);assert stats==cover['statistics']
 print('EXACTLY VERIFIED left-cubic chord identities and',stats['leaves'],'independently checked rational rectangles.',flush=True)
 return stats

def reasons(w):
 ps=patterns(w);sig=ps[0][0]
 if sig[5]==-1 and sig[11]==-1 and sig[0]==-1 and [p[3]for p in ps[:3]]==[1,-1,1]:
  return ['left_cubic_chord']
 return []

if __name__=='__main__':
 check_certificates()
 known=json.loads((ROOT/'data/certificates_n6.json').read_text());unknown=json.loads((ROOT/'data/unresolved_n6.json').read_text())
 assert not any(reasons(w)for w in known)
 print('newly matched:',[w for w in unknown if reasons(w)])
