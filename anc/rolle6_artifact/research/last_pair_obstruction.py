"""Exact last-pair exclusion: matched left-minimum and right-maximum regimes."""
from pathlib import Path
from fractions import Fraction as F
import json
from right_fourth_minimum_obstruction import ONE,U,K,Y,add,scale,mul,power,direct_minimum
from bernstein3 import bernstein,split
from recovery_geometry import check_domain_reduction
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]
def check_certificates(path=None):
 data=json.loads(Path(path or ROOT/'data/last_pair_certificate.json').read_text());parse=lambda row:{tuple(e):F(c)for e,c in row}
 check_domain_reduction();t=K
 D=add(mul(mul(U,t),add(scale(t,3),scale(ONE,-2))),scale(t,2),scale(ONE,-1))
 N=mul(add(scale(mul(U,t),2),ONE),add(mul(U,add(ONE,scale(power(t,2),-2))),scale(t,-2),ONE))
 DN=add(scale(D,2),N);TN=add(mul(add(ONE,U),D),N)
 HN=add(mul(D,add(power(U,2),U,scale(ONE,-1))),mul(N,add(U,scale(ONE,-1))))
 assert all(i>=1 for i,j,k in HN);Hb={(i-1,j,k):c for(i,j,k),c in HN.items()}
 W=add(mul(add(ONE,scale(U,-1)),power(DN,2)),scale(mul(D,DN),F(-1,5)),scale(power(D,2),F(-2,3)))
 ec={6:scale(D,2),5:scale(DN,-3),4:scale(HN,-5),3:scale(mul(U,TN),10)}
 E=lambda x:add(*(mul(c,power(x,j))for j,c in ec.items()))
 EP=lambda x:add(*(scale(mul(c,power(x,j-1)),j)for j,c in ec.items()))
 UT=mul(U,t);UV=mul(U,Y);qa=EP(scale(UT,-1));qx=EP(scale(UV,-1))
 assert all(i>=4 for i,j,k in qa);QA={(i-4,j,k):c for(i,j,k),c in qa.items()}
 assert all(i>=4 for i,j,k in qx);QX={(i-4,j,k):c for(i,j,k),c in qx.items()}
 lhs=add(E(scale(UV,-1)),mul(UV,qa));assert all(i>=5 and k>=1 for i,j,k in lhs)
 Gbar={(i-5,j,k-1):c for(i,j,k),c in lhs.items()}
 Z=add(scale(mul(mul(DN,Y),Gbar),3),scale(mul(QA,Hb),-1))
 lhsq=add(E(scale(UV,-1)),mul(UV,qx));assert all(i>=5 for i,j,k in lhsq)
 Gq={(i-5,j,k):c for(i,j,k),c in lhsq.items()}
 Zq=add(scale(mul(DN,Gq),3),scale(mul(QX,Hb),-1))
 first={'-D':scale(D,-1),'2D-N':add(scale(D,2),scale(N,-1)),'-Hb':scale(Hb,-1),'-W':scale(W,-1),'Z':Z,'qa-qx':add(QA,scale(QX,-1)),'-qx':scale(QX,-1),'Zq':Zq,'Zq_cubic':add(scale(mul(mul(U,power(DN,2)),Zq),9),scale(mul(TN,power(Hb,3)),10))}
 fover=add(*(mul(c,power(Y,j-1))for j,c in ec.items()))
 sec1=add(*(mul(c,add(*(power(Y,k)for k in range(j))))for j,c in ec.items()))
 right={'-D':scale(D,-1),'2D-N':add(scale(D,2),scale(N,-1)),'-Hb':scale(Hb,-1),'W':W,'G0':add(qa,scale(fover,-1)),'G1':add(qa,scale(sec1,-1))}
 assert add(scale(mul(D,power(scale(UT,-1),3)),8),scale(mul(DN,power(UT,2)),-6),scale(mul(HN,UT),4),scale(mul(U,TN),2))=={}
 assert mul(D,add(UT,ONE,U))==add(mul(mul(UT,add(ONE,scale(t,-1))),add(UT,ONE)),scale(N,-1))
 assert F(1,2)*(F(3,4)**2-F(3,4)+1)+2*F(3,4)-2==F(-3,32)
 assert EP(scale(U,F(-3,2)))==scale(mul(power(U,4),add(scale(N,5),mul(D,add(scale(U,14),scale(ONE,10))))),F(-27,16))
 sq=scale(mul(t,power(add(scale(power(Y,2),3),scale(ONE,-5)),2)),F(2,9))
 fl={n:parse(row)for n,row in data['first_local'].items()};L,B,J,C,A=[fl[n]for n in ['L','B','E','C','AA']]
 assert add(Hb,scale(D,-1))==mul(add(U,scale(ONE,-1)),L)
 assert mul(L,add(Gbar,scale(sq,-1)))==add(mul(D,J),mul(Hb,B))
 assert C==add(scale(mul(mul(DN,Y),B),3),scale(mul(L,QA),-1))
 assert A==add(scale(mul(mul(power(DN,2),Y),J),3),mul(add(DN,scale(L,F(-1,5))),C))
 W0=add(mul(add(ONE,scale(U,-1)),DN),scale(D,F(-1,5)))
 assert W==add(mul(DN,W0),scale(power(D,2),F(-2,3)))
 assert mul(mul(DN,L),Z)==add(scale(mul(mul(mul(power(DN,2),Y),L),sq),3),mul(D,A),scale(mul(mul(L,W0),C),-1))
 firstbox=[(F(99,100),F(1)),(F(577,1000),F(578,1000)),(F(128,100),F(130,100))]
 for p in [L,scale(C,-1),A]:
  m=direct_minimum(p,firstbox);assert m>0 and m==min(bernstein(p,firstbox)[1].values())
 rl={n:parse(row)for n,row in data['right_local'].items()};Br,Jr,Ar=[rl[n]for n in ['B','J','A']]
 assert add(right['G0'],scale(sq,-1))==add(mul(D,Ar),mul(add(U,scale(ONE,-1)),Br))
 assert Jr==add(mul(DN,Ar),scale(Br,F(-2,5)))
 W1=add(mul(add(ONE,scale(U,-1)),DN),scale(D,F(-2,5)))
 assert add(W,scale(mul(DN,W1),-1))==scale(mul(D,add(scale(N,3),scale(D,-4))),F(1,15))
 assert mul(DN,add(right['G0'],scale(sq,-1)))==add(mul(D,Jr),scale(mul(W1,Br),-1))
 rightbox=[(F(99,100),F(101,100)),firstbox[1],firstbox[2]]
 for p in [Br,Jr]:
  m=direct_minimum(p,rightbox);assert m>0 and m==min(bernstein(p,rightbox)[1].values())
 def cover(name,polydict,box,localbox):
  c=data[name];assert c['names']==list(polydict)and c['box']==[[str(a),str(b)]for a,b in box]
  polys=list(polydict.values());bs=[bernstein(p,box)for p in polys]
  assert [parse(row)for row in c['polynomials']]==polys
  stats={'nodes':0,'leaves':0,'depth':0,'types':[0]*(len(polys)+1)}
  def walk(tr,bs,box,depth):
   stats['nodes']+=1;stats['depth']=max(stats['depth'],depth)
   if type(tr)is int:
    assert 0<=tr<=len(polys)
    if tr==len(polys):assert all(a>=c and b<=d for(a,b),(c,d)in zip(box,localbox))
    else:
     m=min(bs[tr][1].values());assert m==direct_minimum(polys[tr],box)
     assert m>0 if tr==3 else m>=0
    stats['leaves']+=1;stats['types'][tr]+=1;return
   assert type(tr)is list and len(tr)==3 and type(tr[0])is int and tr[0]in range(3)
   ax=tr[0];mid=sum(box[ax])/2;pairs=[split(p,ax)for p in bs]
   for side in range(2):
    b=box[:];b[ax]=(box[ax][0],mid)if side==0 else(mid,box[ax][1]);walk(tr[side+1],[p[side]for p in pairs],b,depth+1)
  walk(c['tree'],bs,box,0);assert stats==c['statistics'];return stats
 a=cover('first',first,[(F(3,4),F(1)),(F(1,2),F(3,4)),(F(1),F(3,2))],firstbox)
 b=cover('right',right,[(F(3,5),F(2)),(F(1,2),F(3,4)),(F(1),F(2))],rightbox)
 print('EXACTLY VERIFIED last-pair obstruction:',a['leaves'],'left-minimum boxes and',b['leaves'],'right-maximum boxes; both limiting squares and boundaries checked.',flush=True)
 return a,b

def reasons(w):
 p=patterns(w);sig=p[0][0]
 if sig[5]<0 and sig[11]<0 and sig[1]<0 and sig[7]<0 and sig[3]>0 and p[0][3]>0 and p[1][3]<0 and p[2][3]<0:return ['last_pair_split_extrema']
 return []
if __name__=='__main__':check_certificates()
