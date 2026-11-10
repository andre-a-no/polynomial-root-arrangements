"""Middle-interval and convexity obstructions for real-rooted sextics.

Proof: PROOF_MIDDLE_AND_CONVEXITY.md. Arithmetic checks use only the
standard library. They check the algebra in the proof; the real-variable
inequalities and convexity argument are justified in the proof itself.
"""
from fractions import Fraction as F
from analytic_obstruction import add, scale, mul, power, derivative, evaluate_univariate


def check_identities():
    one={(0,0,0):1};u={(1,0,0):1};S={(0,1,0):1};T={(0,0,1):1}
    A={(4,0,0):15,(3,0,0):-30,(2,0,0):5,(1,0,0):10,(0,0,0):-2}
    D={(6,0,0):12,(5,0,0):-22,(4,0,0):4,(3,0,0):2,(2,0,0):5,(1,0,0):-2}
    N=add(scale(mul(T,A),-1),D)
    den=mul(u,add(scale(u,3),scale(one,-2)))
    E60=add(scale(one,-1),scale(S,2),scale(T,-5),scale(power(u,5),-12),
        scale(mul(add(one,S),power(u,4)),15),
        scale(mul(add(S,T),power(u,3)),-20),scale(mul(T,power(u,2)),30))
    # Compute 60*(integral G - G(u)) independently from G'=R.
    rcoeff=[{},scale(T,-1),add(S,T),scale(add(one,S),-1),one]
    actual={}
    for k,coeff in enumerate(rcoeff):
        actual=add(actual,scale(coeff,F(1,(k+1)*(k+2))),
            scale(mul(coeff,power(u,k+1)),F(-1,k+1)))
    assert scale(actual,60)==E60
    Rprime=add(scale(power(u,3),4),scale(mul(add(one,S),power(u,2)),-3),
        scale(mul(add(S,T),u),2),scale(T,-1))
    odd=add(scale(u,2),scale(one,-1))
    M=add(scale(power(u,4),-15),scale(power(u,3),20),scale(one,-2))
    assert not add(mul(den,E60),scale(mul(odd,N),-1),scale(mul(M,Rprime),-1))
    four_mu_diff=add(scale(u,4),scale(one,-1),scale(S,-1))
    H=add(scale(T,-1),scale(power(u,2),4),scale(u,-2))
    assert not add(mul(den,four_mu_diff),scale(mul(odd,H),-1),scale(Rprime,-1))
    assert derivative(A)==scale(mul(odd,
        add(scale(power(u,2),3),scale(u,-3),scale(one,-1))),10)
    assert evaluate_univariate(A,F(2,3))==F(26,27)
    v=add(u,scale(one,F(-1,2)))
    expanded_D=add(scale(power(v,6),12),scale(power(v,5),14),
        scale(power(v,4),-6),scale(power(v,3),-15),
        scale(power(v,2),F(-9,4)),scale(v,F(15,8)),scale(one,F(1,4)))
    assert D==expanded_D
    assert F(15,8)-F(9,4)*F(1,6)-15*F(1,36)-6*F(1,216)==F(19,18)
    # A cubic with roots 0,1,S has derivative -1/4 at 1/2,
    # independently of its third root S. This checks both midpoint claims
    # after separately normalizing each of its adjacent root intervals.
    cubic=mul(mul(u,add(u,scale(one,-1))),add(u,scale(S,-1)))
    at_half={}
    for m,c in derivative(cubic).items():
        at_half=add(at_half,{(0,m[1],m[2]):c*F(1,2)**m[0]})
    assert at_half==scale(one,F(-1,4))
    print('EXACTLY VERIFIED middle-interval identities, rational positivity bounds and cubic midpoint identity',flush=True)


def reasons(w):
    assert len(w)==21 and all(w.count(str(i))==6-i for i in range(6))
    pos=[[j for j,c in enumerate(w)if c==str(i)]for i in range(6)]
    sign=lambda row,x:(-1)**sum(j>x for j in pos[row])
    out=[]
    if pos[3][1]<pos[5][0]:
        if sign(0,pos[2][1])<0<sign(0,pos[2][2]) and sign(1,pos[3][1])<0:
            out.append('middle_left_of_mean')
    else:
        if sign(0,pos[2][1])>0>sign(0,pos[2][2]) and sign(1,pos[3][1])>0:
            out.append('middle_right_of_mean')
    if (sign(1,pos[4][0])>0 and sign(2,pos[4][0])>0
            and sign(0,pos[3][0])>0>sign(0,pos[3][1])):
        out.append('convex_first')
    if (sign(1,pos[4][1])<0 and sign(2,pos[4][1])>0
            and sign(0,pos[3][1])<0<sign(0,pos[3][2])):
        out.append('concave_last')
    return out


if __name__=='__main__':
    check_identities()
