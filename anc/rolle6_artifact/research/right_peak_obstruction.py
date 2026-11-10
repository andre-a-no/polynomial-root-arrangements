"""Exact obstruction forcing the fourth critical maximum to be negative."""
from pathlib import Path
from fractions import Fraction as F
import json
from right_fourth_minimum_obstruction import ONE,U,K,Y,add,scale,mul,power,translate,direct_minimum
from bernstein3 import bernstein,split
from inflection_obstruction import patterns
ROOT=Path(__file__).resolve().parents[1]
NAMES=['-N','N-2D','-H','x-T','mean','right4','left']
def check_certificates(path=None):
 d=json.loads(Path(path or ROOT/'data/right_peak_certificate.json').read_text())
 src={n:{tuple(e):F(c)for e,c in d['polynomials'][n]}for n in NAMES}
 r=K;T=add(ONE,U,scale(r,2));s=add(ONE,r);H=add(power(U,2),U,scale(mul(U,r),2),scale(r,-2),scale(ONE,-1))
 cf={6:scale(ONE,F(1,30)),5:scale(s,F(-1,10)),4:scale(H,F(-1,12)),3:scale(mul(U,T),F(1,6))}
 X=add(ONE,mul(add(U,scale(r,2)),Y))
 tail=scale(add(*(scale(mul(c,power(U,j-1)),-(-1)**j)for j,c in cf.items()),*(scale(mul(c,power(X,j-1)),-1)for j,c in cf.items())),60)
 tail=translate(tail,(F(2),F(0),F(0)));tb=[]
 for j in range(6):
  p={(0,k,l):c for(i,k,l),c in tail.items()if i==j};box=[(F(0),F(1))]*3
  m=direct_minimum(p,box);assert m>0 and m==min(bernstein(p,box)[1].values());tb.append(str(m))
 assert tb==['1183/5','2734/5','474','196','39','3']
 k=K;D=add(U,scale(ONE,2),scale(k,3));N=add(scale(power(k,2),6),scale(k,6),scale(power(U,2),-1),scale(U,-1),ONE)
 DN=add(scale(D,2),N);TN=add(mul(add(ONE,U),D),N)
 HN=add(mul(D,add(power(U,2),U,scale(ONE,-1))),mul(N,add(U,scale(ONE,-1))))
 assert src['-N']==scale(N,-1)and src['N-2D']==add(N,scale(D,-2))and src['-H']==scale(HN,-1)
 assert src['x-T']==add(mul(D,Y),scale(TN,-1))
 V=add(ONE,k)
 assert add(scale(mul(D,power(V,2)),12),scale(mul(DN,V),-6),scale(HN,-2))=={}
 ec={6:scale(D,2),5:scale(DN,-3),4:scale(HN,-5),3:scale(mul(U,TN),10)}
 Fx=add(*(mul(c,power(Y,j-1))for j,c in ec.items()))
 Qv=add(*(scale(mul(c,power(V,j-1)),j)for j,c in ec.items()))
 assert src['right4']==add(Qv,scale(Fx,-1))
 assert src['left']==add(*(scale(mul(c,power(U,j-1)),-(-1)**j)for j,c in ec.items()),scale(Fx,-1))
 fourD=scale(D,4)
 Qmean_num=add(*(scale(mul(mul(c,power(DN,j-1)),power(fourD,6-j)),j)for j,c in ec.items()))
 assert scale(mul(D,src['mean']),32)==add(Qmean_num,scale(mul(power(D,5),Fx),-1024))
 # H>0 implies u>3/5, since H(3/5,r)=-1/25-4r/5.
 assert F(3,5)**2+F(3,5)-1==F(-1,25)
 # On u<=2,r<=1: R''(5/2)>=1, hence k<3/2; also T<=5.
 assert 47-26-2*4-2*2-4*2==1
 box=[(F(3,5),F(2)),(F(0),F(3,2)),(F(1),F(5))]
 cover=d['cover'];assert cover['names']==NAMES and cover['box']==[[str(a),str(b)]for a,b in box]
 polys=list(src.values());bs=[bernstein(p,box)for p in polys];stat={'nodes':0,'leaves':0,'depth':0,'types':[0]*7}
 def walk(tr,bs,box,depth):
  stat['nodes']+=1;stat['depth']=max(stat['depth'],depth)
  if type(tr)is int:
   assert 0<=tr<7;m=min(bs[tr][1].values());assert m>=0 and m==direct_minimum(polys[tr],box)
   stat['leaves']+=1;stat['types'][tr]+=1;return
  assert type(tr)is list and len(tr)==3 and type(tr[0])is int and tr[0]in range(3)
  ax=tr[0];pairs=[split(p,ax)for p in bs];mid=sum(box[ax])/2
  for side in range(2):
   b=box[:];b[ax]=(box[ax][0],mid)if side==0 else(mid,box[ax][1]);walk(tr[side+1],[p[side]for p in pairs],b,depth+1)
 walk(cover['tree'],bs,box,0);assert stat==cover['statistics']
 print('EXACTLY VERIFIED right-peak obstruction:',stat['leaves'],'independent rational boxes and six positive tail coefficients.',flush=True)
 return stat

def reasons(w):
 p=patterns(w);sig=p[0][0]
 if sig[5]<0 and sig[11]>0 and sig[10]<0 and sig[13]<0 and p[0][3]>0 and p[1][3]<0:return ['negative_fourth_critical_maximum']
 return []
if __name__=='__main__':check_certificates()
