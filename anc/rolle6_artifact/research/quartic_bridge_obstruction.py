"""Middle critical root of a quartic: between mean and inner-root midpoint."""
import json
from metric_obstruction import ROOT,constraint_forms as base_forms
from subfactor_sign_rules import propagate
from root_subset_coordinates import FULL,SUBSETS,D,coordinates,linear,verify_sparse

def constraint_forms(w):
    E,B=base_forms(w);Q=coordinates(w);info=propagate(w);assert 'conflict'not in info;center=w.index('5')
    for S in [FULL]+SUBSETS:
        n=len(S)
        if n<4:continue
        j=n-4;roots=Q[S][j]
        if S==FULL:
            mid=[i for i,c in enumerate(w)if c=='3'][1];sign=1 if mid>center else -1
        else:
            a,b=info['bounds'][S,j+1][1];c,d=info['bounds'][S,n-1][0]
            sign=1 if a>=d else (-1 if b<=c else 0)
        if sign:
            B.append(linear([(sign,roots[1]),(-sign,roots[0]),(-sign,roots[3]),(sign,roots[2])]))
            B.append(linear([(sign,roots[1]),(sign,roots[2]),(-2*sign,Q[S][j+1][1])]))
    return E,B

def verify_certificate(w,record):verify_sparse(*constraint_forms(w),record)


def check_identities():
    from analytic_obstruction import add,mul,scale
    a={(1,0,0):1};b={(0,1,0):1};c={(0,0,1):1};d=scale(add(a,b,c),-1)
    roots=[a,b,c,d]
    e3=add(*(mul(mul(roots[i],roots[j]),roots[k])for i in range(4)for j in range(i+1,4)for k in range(j+1,4)))
    assert not add(scale(e3,-1),scale(mul(mul(add(a,b),add(a,c)),add(b,c)),-1))
    from fractions import Fraction
    h=scale(add(b,c),Fraction(1,2))
    at_h=add(*(mul(mul(add(h,scale(roots[i],-1)),add(h,scale(roots[j],-1))),add(h,scale(roots[k],-1)))for i in range(4)for j in range(i+1,4)for k in range(j+1,4)))
    rhs=mul(mul(add(scale(h,2),scale(a,-1),scale(d,-1)),add(h,scale(b,-1))),add(h,scale(c,-1)))
    assert not add(at_h,scale(rhs,-1))
    print('EXACTLY VERIFIED two quartic midpoint identities',flush=True)

def check_certificates():
    check_identities()
    data=json.loads((ROOT/'data/quartic_bridge_rejected_n6.json').read_text())
    for w,record in data.items():verify_certificate(w,record)
    result=set(data);assert {w[::-1]for w in result}==result
    print('EXACTLY VERIFIED',len(result),'quartic middle-root certificates',flush=True)
    return result

if __name__=='__main__':check_certificates()
