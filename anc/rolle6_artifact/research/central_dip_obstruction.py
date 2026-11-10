"""Exact conflict between the first minimum and a central negative value."""
from pathlib import Path
from fractions import Fraction as F
import json
from right_fourth_minimum_obstruction import ONE,U,K,Y,add,scale,mul,power,direct_minimum
from bernstein3 import bernstein
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]
BOX=[(F(2,5),F(1)),(F(1,2),F(3,4)),(F(1),F(3,2))]
def geometry():
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
 lhs=add(E(scale(UV,-1)),mul(UV,qx));assert all(i>=5 for i,j,k in lhs)
 Gq={(i-5,j,k):c for(i,j,k),c in lhs.items()}
 base=[scale(D,-1),scale(N,-1),add(N,scale(D,-2)),Hb,scale(add(HN,scale(N,3)),-1),add(QA,scale(QX,-1)),scale(QX,-1),Gq]
 cubic=lambda c:add(scale(mul(D,power(scale(ONE,c),3)),8),scale(mul(DN,power(scale(ONE,c),2)),-6),scale(mul(HN,scale(ONE,c)),-4),scale(mul(U,TN),2))
 assert cubic_poly_at_a(D,DN,HN,U,TN,UT)=={}
 assert mul(D,add(UT,ONE,U))==add(mul(mul(UT,add(ONE,scale(t,-1))),add(UT,ONE)),scale(N,-1))
 assert EP(scale(U,F(-3,2)))==scale(mul(power(U,4),add(scale(N,5),mul(D,add(scale(U,14),scale(ONE,10))))),F(-27,16))
 central=lambda c:add(scale(E(scale(ONE,c)),1/c),scale(qx,-1))
 second=lambda w:add(scale(mul(D,power(scale(ONE,w),2)),24),scale(mul(DN,scale(ONE,w)),-12),scale(HN,-4))
 upper=lambda w:add(qx,scale(EP(scale(ONE,w)),-1))
 return base,cubic,central,second,upper
def cubic_poly_at_a(D,DN,HN,U,TN,UT):
 return add(scale(mul(D,power(scale(UT,-1),3)),8),scale(mul(DN,power(UT,2)),-6),scale(mul(HN,UT),4),scale(mul(U,TN),2))
BASE,CUBIC,CENTRAL,SECOND,UPPER=geometry()
def check_domain_reduction():
 r=K;T=add(ONE,U,scale(r,2));H=add(power(U,2),U,scale(mul(U,r),2),scale(r,-2),scale(ONE,-1))
 cf={6:scale(ONE,F(1,30)),5:scale(add(ONE,r),F(-1,10)),4:scale(H,F(-1,12)),3:scale(mul(U,T),F(1,6))}
 E=lambda x:add(*(mul(c,power(x,j))for j,c in cf.items()))
 Q=lambda x:add(*(scale(mul(c,power(x,j-1)),j)for j,c in cf.items()))
 p=add(scale(E(scale(ONE,F(1,2))),2),scale(Q(scale(U,-1)),-1))
 box=[(F(0),F(2,5)),(F(0),F(1)),(F(0),F(1))]
 assert direct_minimum(p,box)==min(bernstein(p,box)[1].values())==F(2143,1000000)
 assert Q(scale(U,F(-3,2)))==scale(mul(power(U,4),add(scale(r,5),scale(U,7),scale(ONE,5))),F(-9,160))
def check_certificates(path=None):
 check_domain_reduction();d=json.loads(Path(path or ROOT/'data/central_dip_certificate.json').read_text())
 assert d['box']==[[str(a),str(b)]for a,b in BOX]
 stats={'nodes':0,'leaves':0,'depth':0,'types':[0]*10};cache=set()
 def walk(tr,box,depth):
  stats['nodes']+=1;stats['depth']=max(stats['depth'],depth)
  if isinstance(tr,dict):
   i=tr['leaf'];assert type(i)is int and 0<=i<10
   if i<8:p=BASE[i]
   elif i==8:
    c=F(tr['c_lower']);assert F(1,2)<=c<F(1)
    if c!=F(1,2):
     key=('c',tuple(box[:2]),c)
     if key not in cache:
      pbox=box[:2]+[(F(0),F(1))];m=direct_minimum(CUBIC(c),pbox)
      assert m>=0 and m==min(bernstein(CUBIC(c),pbox)[1].values());cache.add(key)
    p=CENTRAL(c)
   else:
    w=F(tr['vplus_lower']);assert F(1)<=w<F(2)
    if w!=F(1):
     key=('vplus',tuple(box[:2]),w)
     if key not in cache:
      pbox=box[:2]+[(F(0),F(1))];p=scale(SECOND(w),-1);m=direct_minimum(p,pbox)
      assert m>=0 and m==min(bernstein(p,pbox)[1].values());cache.add(key)
    p=UPPER(w)
   m=direct_minimum(p,box);assert m>=0 and m==min(bernstein(p,box)[1].values()),(box,i,m)
   stats['leaves']+=1;stats['types'][i]+=1;return
  assert type(tr)is list and len(tr)==3 and type(tr[0])is int and tr[0]in range(3)
  ax=tr[0];mid=sum(box[ax])/2
  for side in range(2):
   b=box[:];b[ax]=(box[ax][0],mid)if side==0 else(mid,box[ax][1]);walk(tr[side+1],b,depth+1)
 walk(d['tree'],BOX,0);assert stats==d['statistics']
 print('EXACTLY VERIFIED central-dip obstruction:',stats['leaves'],'rational boxes and certified central-root bounds.',flush=True)
 return stats
def reasons(w):
 p=patterns(w);s=p[0][0]
 if s[5]>0 and s[11]>0 and s[14]<0 and s[1]<0 and s[6]<0 and s[13]>0 and p[1][3]>0:return ['central_dip_first_minimum']
 return []
if __name__=='__main__':check_certificates()
