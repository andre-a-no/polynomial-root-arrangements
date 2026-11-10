"""A negative left-cubic value prevents a positive fourth critical maximum."""
from pathlib import Path
from fractions import Fraction as F
import json
from right_fourth_minimum_obstruction import ONE,U,K,Y,add,scale,mul,power,direct_minimum
from bernstein3 import bernstein
from inflection_obstruction import patterns
from recovery_geometry import HH,COEFF,f_at,q_at,cubic_at,secant,check_domain_reduction
ROOT=Path(__file__).resolve().parents[1]
BASE=[scale(HH,-1)]
def chord(x,z):return add(*(mul(co,add(*(mul(power(x,j-1-k),power(z,k))for k in range(j))))for j,co in COEFF.items()))
def polynomials(anchors):
 a=anchors[0];A=scale(ONE,a);qa=q_at(A);left=chord(A,scale(U,-1))
 return BASE+[add(qa,scale(q_at(ONE),-1))]+[add(q,scale(secant(z),-1))for q in [qa,left]for z in [{},ONE,A]]
def anchor_checks(box,anchors):
 a=anchors[0];assert a<0
 assert a==-box[0][0]/2 or direct_minimum(cubic_at(a),[box[0],box[1],(F(0),F(1))])>=0

def check_certificates(path=None):
 check_domain_reduction()
 d=json.loads(Path(path or ROOT/'data/left_cubic_negative_certificate.json').read_text());box=[(F(3,5),F(2)),(F(0),F(1,2)),(F(1),F(2))]
 assert d['box']==[[str(a),str(b)]for a,b in box]
 stats={'nodes':0,'leaves':0,'depth':0,'types':[0]*8};cache={}
 def walk(tr,box,depth):
  stats['nodes']+=1;stats['depth']=max(stats['depth'],depth)
  if isinstance(tr,dict):
   i=tr['leaf'];assert type(i)is int and 0<=i<8
   if i==0:p=BASE[0]
   else:
    anchors=tuple(F(x)for x in tr['anchors']);key=(tuple(box[:2]),anchors)
    if key not in cache:anchor_checks(box,anchors);cache[key]=True
    p=polynomials(anchors)[i]
   m=direct_minimum(p,box);assert m>=0 and m==min(bernstein(p,box)[1].values())
   stats['leaves']+=1;stats['types'][i]+=1;return
  assert type(tr)is list and len(tr)==3 and type(tr[0])is int and tr[0]in range(3)
  ax=tr[0];mid=sum(box[ax])/2
  for side in range(2):
   b=box[:];b[ax]=(box[ax][0],mid)if side==0 else(mid,box[ax][1]);walk(tr[side+1],b,depth+1)
 walk(d['tree'],box,0);assert stats==d['statistics']
 print('EXACTLY VERIFIED left-cubic negative-value obstruction:',stats['leaves'],'rational boxes.',flush=True)
 return stats

def reasons(w):
 p=patterns(w);sig=p[0][0]
 if sig[5]<0 and sig[11]<0 and sig[1]<0 and sig[0]<0 and p[0][3]>0 and p[1][3]<0 and p[2][3]<0:return ['left_cubic_negative_fourth_maximum']
 return []
if __name__=='__main__':check_certificates()
