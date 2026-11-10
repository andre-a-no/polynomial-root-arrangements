"""Exact obstruction: positive rebound at the mean forces a positive first minimum."""
from pathlib import Path
from fractions import Fraction as F
import json
from right_fourth_minimum_obstruction import ONE,U,K,Y,add,scale,mul,power,translate,direct_minimum
from bernstein3 import bernstein,split
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1];NAMES=['H','-S','-C','G']
def check_certificates(path=None):
 d=json.loads(Path(path or ROOT/'data/mean_rebound_certificate.json').read_text())
 src={n:{tuple(e):F(c)for e,c in d['polynomials'][n]}for n in NAMES}
 r=K;s=add(ONE,r);T=add(ONE,U,scale(r,2));H=add(power(U,2),U,scale(mul(U,r),2),scale(r,-2),scale(ONE,-1));mu=scale(s,F(1,2))
 cf={6:scale(ONE,F(1,30)),5:scale(s,F(-1,10)),4:scale(H,F(-1,12)),3:scale(mul(U,T),F(1,6))}
 ev=lambda x:add(*(mul(c,power(x,i))for i,c in cf.items()))
 qv=lambda x:add(*(scale(mul(c,power(x,i-1)),i)for i,c in cf.items()))
 f1=ev(ONE);yy=add(U,Y)
 assert src['H']==H
 assert scale(add(ev(mu),scale(mul(mu,f1),-1)),960)==scale(mul(mul(mu,add(mu,scale(ONE,-1))),src['-S']),-1)
 assert scale(add(qv(scale(U,-1)),scale(f1,-1)),60)==scale(src['-C'],-1)
 assert scale(add(ev(scale(yy,-1)),mul(yy,f1)),60)==mul(mul(yy,add(yy,ONE)),src['G'])
 # r>=3 is impossible if S>0. Every coefficient in r-3 is positive.
 tail=translate(src['-S'],(F(0),F(3),F(0)));rb=[]
 for j in range(5):
  p={(i,0,k):c for(i,l,k),c in tail.items()if l==j};b=[(F(0),F(1))]*3
  m=direct_minimum(p,b);assert m>0 and m==min(bernstein(p,b)[1].values());rb.append(str(m))
 assert rb==['848','1400','700','140','10']
 # On u<=3/5,1<=r<=3, C>0 is impossible.
 small=[(F(0),F(3,5)),(F(1),F(3)),(F(0),F(1))]
 m=direct_minimum(src['-C'],small);assert m==F(15031,3125)==min(bernstein(src['-C'],small)[1].values())
 # Cover the entire unbounded spatial tail z>=2.
 tail=translate(src['G'],(F(0),F(0),F(2)));zb=[]
 box=[(F(3,5),F(1)),(F(1),F(3)),(F(0),F(2))]
 for j in range(5):
  p={(i,k,0):c for(i,k,l),c in tail.items()if l==j};b=[box[0],box[1],(F(0),F(1))]
  m=direct_minimum(p,b);assert m>0 and m==min(bernstein(p,b)[1].values());zb.append(str(m))
 assert zb==['126667/625','37181/125','3833/25','154/5','2']
 cover=d['cover'];assert cover['names']==NAMES and cover['box']==[[str(a),str(b)]for a,b in box]
 polys=list(src.values());bs=[bernstein(p,box)for p in polys];stat={'nodes':0,'leaves':0,'depth':0,'types':[0]*4}
 def walk(tr,bs,box,depth):
  stat['nodes']+=1;stat['depth']=max(stat['depth'],depth)
  if type(tr)is int:
   assert 0<=tr<4;m=min(bs[tr][1].values());assert m>=0 and m==direct_minimum(polys[tr],box)
   stat['leaves']+=1;stat['types'][tr]+=1;return
  assert type(tr)is list and len(tr)==3 and type(tr[0])is int and tr[0]in range(3)
  ax=tr[0];pairs=[split(p,ax)for p in bs];mid=sum(box[ax])/2
  for side in range(2):
   b=box[:];b[ax]=(box[ax][0],mid)if side==0 else(mid,box[ax][1]);walk(tr[side+1],[p[side]for p in pairs],b,depth+1)
 walk(cover['tree'],bs,box,0);assert stat==cover['statistics']
 print('EXACTLY VERIFIED mean-rebound obstruction:',stat['leaves'],'rational boxes; small-u exclusion, parameter tail, and complete spatial tail.',flush=True)
 return stat

def reasons(w):
 p=patterns(w);sig=p[0][0]
 if sig[5]>0 and sig[11]<0 and sig[9]>0 and p[1][3]>0 and p[2][3]<0:return ['mean_rebound_first_minimum']
 return []
if __name__=='__main__':check_certificates()
