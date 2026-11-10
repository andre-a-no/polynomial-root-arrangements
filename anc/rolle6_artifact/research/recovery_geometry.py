"""Shared exact geometry for the four terminal words, recovered 2026-11-10."""
from fractions import Fraction as F
from right_fourth_minimum_obstruction import ONE,U,K,Y,add,scale,mul,power,translate,direct_minimum
h=K
HH=add(mul(h,add(power(U,2),scale(U,-1),ONE)),scale(U,2),scale(ONE,-2))
TT=add(scale(ONE,2),mul(h,add(U,scale(ONE,-1))))
COEFF={6:scale(h,F(1,30)),5:scale(ONE,F(-1,10)),4:scale(HH,F(-1,12)),3:scale(mul(U,TT),F(1,6))}
def f_at(x):return add(*(mul(co,power(x,j))for j,co in COEFF.items()))
def q_at(x):return add(*(scale(mul(co,power(x,j-1)),j)for j,co in COEFF.items()))
def cubic_at(c):
 x=scale(ONE,c)
 return add(scale(mul(h,power(x,3)),4),scale(power(x,2),-6),scale(mul(HH,x),-2),mul(U,TT))
def secant(z):return add(*(mul(co,add(*(mul(power(Y,j-1-k),power(z,k))for k in range(j))))for j,co in COEFF.items()))
def check_domain_reduction():
 r=K;T=add(ONE,U,scale(r,2));H=add(power(U,2),U,scale(mul(U,r),2),scale(r,-2),scale(ONE,-1))
 cf={5:scale(ONE,F(1,5)),4:scale(add(ONE,r),F(-1,2)),3:scale(H,F(-1,3)),2:scale(mul(U,T),F(1,2))}
 Q=lambda x:add(*(mul(c,power(x,j))for j,c in cf.items()))
 tail=translate(scale(add(Q(scale(U,F(-1,2))),scale(Q(ONE),-1)),480),(F(2),F(1),F(0)))
 assert all(c>0 for c in tail.values())and tail[(0,0,0)]==2048
 assert scale(Q(scale(ONE,2)),15)==scale(add(scale(power(U,2),10),scale(U,10),scale(mul(U,r),20),scale(r,40),scale(ONE,-16)),-1)
