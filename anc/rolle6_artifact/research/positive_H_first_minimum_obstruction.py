"""Exact exclusion from a positive second inflection and the first cubic slope."""
from pathlib import Path
from fractions import Fraction as F
import json
from right_fourth_minimum_obstruction import ONE,U,K,Y,add,scale,mul,power,translate,direct_minimum
from bernstein3 import bernstein,split
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]
BOX=[(F(3,5),F(2)),(F(1,2),F(3,4)),(F(1),F(3,2))]
LOCAL_BOX=[(F(99,100),F(101,100)),(F(577,1000),F(578,1000)),(F(128,100),F(130,100))]
def polynomials():
 t=K
 D=add(mul(mul(U,t),add(scale(t,3),scale(ONE,-2))),scale(t,2),scale(ONE,-1))
 N=mul(add(scale(mul(U,t),2),ONE),add(mul(U,add(ONE,scale(power(t,2),-2))),scale(t,-2),ONE))
 DN=add(scale(D,2),N);TN=add(mul(add(ONE,U),D),N)
 HN=add(mul(D,add(power(U,2),U,scale(ONE,-1))),mul(N,add(U,scale(ONE,-1))))
 assert all(i>=1 for i,j,k in HN);Hb={(i-1,j,k):c for(i,j,k),c in HN.items()}
 ec={6:scale(D,2),5:scale(DN,-3),4:scale(HN,-5),3:scale(mul(U,TN),10)}
 E=lambda x:add(*(mul(c,power(x,j))for j,c in ec.items()))
 EP=lambda x:add(*(scale(mul(c,power(x,j-1)),j)for j,c in ec.items()))
 UT=mul(U,t);UV=mul(U,Y);qa=EP(scale(UT,-1));qx=EP(scale(UV,-1))
 assert all(i>=4 for i,j,k in qa);QA={(i-4,j,k):c for(i,j,k),c in qa.items()}
 assert all(i>=4 for i,j,k in qx);QX={(i-4,j,k):c for(i,j,k),c in qx.items()}
 lhs=add(E(scale(UV,-1)),mul(UV,qa));assert all(i>=5 and k>=1 for i,j,k in lhs)
 G={(i-5,j,k-1):c for(i,j,k),c in lhs.items()}
 # R'(-ut)=0 and the exact denominator identity from logarithmic differentiation.
 assert add(scale(mul(D,power(scale(UT,-1),3)),8),scale(mul(DN,power(UT,2)),-6),scale(mul(HN,UT),4),scale(mul(U,TN),2))=={}
 assert mul(D,add(UT,ONE,U))==add(mul(mul(UT,add(ONE,scale(t,-1))),add(UT,ONE)),scale(N,-1))
 assert EP(scale(U,F(-3,2)))==scale(mul(power(U,4),add(scale(N,5),mul(D,add(scale(U,14),scale(ONE,10))))),F(-27,16))
 return {'-D':scale(D,-1),'-N':scale(N,-1),'-Hb':scale(Hb,-1),'G':G,'qa-qx':add(QA,scale(QX,-1)),'-qx':scale(QX,-1)},D,Hb,G

def check_domain_reduction():
 r=K;T=add(ONE,U,scale(r,2));H=add(power(U,2),U,scale(mul(U,r),2),scale(r,-2),scale(ONE,-1))
 cf={5:scale(ONE,F(1,5)),4:scale(add(ONE,r),F(-1,2)),3:scale(H,F(-1,3)),2:scale(mul(U,T),F(1,2))}
 Q=lambda x:add(*(mul(c,power(x,j))for j,c in cf.items()))
 tail=translate(scale(add(Q(scale(U,F(-1,2))),scale(Q(ONE),-1)),480),(F(2),F(0),F(0)))
 assert all(c>0 for c in tail.values())and tail[(0,0,0)]==1408
 assert Q(scale(U,F(-3,2)))==scale(mul(power(U,4),add(scale(r,5),scale(U,7),scale(ONE,5))),F(-9,160))

def check_certificates(path=None):
 check_domain_reduction();d=json.loads(Path(path or ROOT/'data/positive_H_first_minimum_certificate.json').read_text())
 parse=lambda row:{tuple(e):F(c)for e,c in row}
 pd,D,Hb,G=polynomials();ps=list(pd.values())
 assert d['box']==[[str(a),str(b)]for a,b in BOX]and d['names']==list(pd)
 assert [parse(row)for row in d['polynomials']]==ps
 L,B,E=[parse(d['local'][n])for n in ['L','B','E']]
 sq=scale(mul(K,power(add(scale(power(Y,2),3),scale(ONE,-5)),2)),F(2,9))
 assert add(Hb,scale(D,-1))==mul(add(U,scale(ONE,-1)),L)
 assert mul(L,add(G,scale(sq,-1)))==add(mul(D,E),mul(Hb,B))
 for p in [L,B,E]:
  m=direct_minimum(p,LOCAL_BOX);assert m>0 and m==min(bernstein(p,LOCAL_BOX)[1].values())
 stats={'nodes':0,'leaves':0,'depth':0,'types':[0]*(len(ps)+1)}
 def walk(tr,bs,box,depth):
  stats['nodes']+=1;stats['depth']=max(stats['depth'],depth)
  if type(tr)is int:
   assert 0<=tr<=len(ps)
   if tr==len(ps):assert all(a>=c and b<=dd for(a,b),(c,dd)in zip(box,LOCAL_BOX))
   else:
    m=min(bs[tr][1].values());assert m>=0 and m==direct_minimum(ps[tr],box)
   stats['leaves']+=1;stats['types'][tr]+=1;return
  assert type(tr)is list and len(tr)==3 and type(tr[0])is int and tr[0]in range(3)
  ax=tr[0];mid=sum(box[ax])/2;pairs=[split(p,ax)for p in bs]
  for side in range(2):
   b=box[:];b[ax]=(box[ax][0],mid)if side==0 else(mid,box[ax][1]);walk(tr[side+1],[p[side]for p in pairs],b,depth+1)
 walk(d['tree'],[bernstein(p,BOX)for p in ps],BOX,0);assert stats==d['statistics']
 print('EXACTLY VERIFIED positive-H first-minimum obstruction:',stats['leaves'],'rational boxes; limiting square and unbounded r included.',flush=True)
 return stats

def reasons(w):
 p=patterns(w);s=p[0][0]
 if s[5]<0 and s[1]<0 and p[1][3]>0:return ['positive_H_first_minimum']
 return []
if __name__=='__main__':check_certificates()
