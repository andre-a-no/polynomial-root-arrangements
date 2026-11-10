"""Resolve simultaneous k=1 / stationary-point boundary degenerations."""
import json
from pathlib import Path
from fractions import Fraction as F
from full_stationary_obstruction import raw_inequalities,witness_polynomials,witness_positive,patterns,pat_key,BOX,dense
from full_stationary_obstruction import check_combination_leaf
from two_parameter_obstruction import add,scale,mul,power,dense_to_sparse,bernstein2,split2
ROOT=Path(__file__).resolve().parents[1]

def chart_data(name):
    one={(0,0):1};K={(1,0):1};T={(0,1):1}
    if name.startswith('bounded_'):
        _,side,bound=name.split('_');assert side in('left','right')and bound in('4','5half')
        D=add(power(K,2),scale(K,-1),one)
        AN=add(scale(one,2),scale(K,-2),scale(power(K,2),-1))
        CN=add(scale(one,-1),scale(K,-2),scale(power(K,2),2));DN=scale(add(AN,CN),-1)
        num,den=(4,1)if bound=='4'else(5,2)
        delta=add(scale(D,num),scale(AN,den))if side=='left'else add(scale(D,num),scale(DN,-den))
        N=mul(delta,T);H=add(scale(D,den),N);return N,H
    base=name.removeprefix('reverse_')
    if base.startswith('right_'):
        q=add(one,scale(K,2),scale(power(K,2),-2));base=base.removeprefix('right_')
    else:q=add(K,scale(one,-1))
    if base=='slope':N=scale(mul(q,T),2);H=add(one,N)
    elif base=='inverse':N=q;H=add(q,T)
    else:raise AssertionError(name)
    if name.startswith('reverse_'):N=add(H,scale(N,-1))
    return N,H

def transform(p,name):
    N,H=chart_data(name);m=max((j for i,j in p),default=0)
    return add(*(scale(mul({(i,0):1},mul(power(N,j),power(H,m-j))),v)for(i,j),v in p.items()))

def check_certificates(path=None):
    data=json.loads(Path(path or ROOT/'data/blowup_certificates.json').read_text())
    assert data['rectangle']==[[1,1],[11,8],[0,1],[1,1]]
    seen=set();leaves=0
    for cert in data['certificates']:
        pat=pat_key(cert);assert pat not in seen;seen.add(pat)
        names=[c['chart']for c in cert['charts']]
        rectangles=[c.get('rectangle',data['rectangle'])for c in cert['charts']]
        full=data['rectangle'];left=[[1,1],[5,4],[0,1],[1,1]];right=[[5,4],[11,8],[0,1],[1,1]]
        allowed_two=[['slope','inverse'],['reverse_slope','reverse_inverse'],
                     ['right_slope','right_inverse'],['reverse_right_slope','reverse_right_inverse']]
        if len(names)==1:
            # P' has five distinct real roots, mean 0 and sum of squares 20.
            # Thus all have |u|<4; its second and fourth roots have |u|<sqrt(6)<5/2.
            assert pat[1]in(0,6)
            side='left'if pat[1]==0 else'right';bound='5half'if pat[2]>0 else'4'
            assert names==['bounded_'+side+'_'+bound]and rectangles==[full]
        elif len(names)==2:assert names in allowed_two and rectangles==[full,full]
        else:
            assert len(names)==4 and names[:2]in allowed_two[:2]and names[2:]in allowed_two[2:]
            assert rectangles==[left,left,right,right]
        # Their open images cover 0<s<1, with overlap:
        # q/(q+1) < 2q/(1+2q), q=k-1>0. Both denominators >0.
        ineq,outside=raw_inequalities(pat)
        for chart in cert['charts']:
            rect=tuple(F(*v)for v in chart.get('rectangle',data['rectangle']))
            initial=[]
            for r in chart['witnesses']:
                raw=[outside]if r['pair']==['outside']else witness_polynomials(ineq,*r['pair'])
                assert len(raw)==len(r['polynomials'])
                for p,pr in zip(raw,r['polynomials']):
                    expected=transform(p,chart['chart']);q=dense_to_sparse(pr['coefficients']);factor={(0,0):1}
                    assert all(type(v)is int for row in pr['coefficients']for v in row)
                    for fr in pr['factors']:
                        fp=dense_to_sparse(fr['coefficients']);n=fr['exponent'];assert type(n)is int and n>0
                        if fr['kind']=='physical':assert fp==scale(outside,-1)
                        elif fr['kind']=='chart_denominator':assert fp==chart_data(chart['chart'])[1]
                        else:
                            assert fr['kind']=='bernstein';bs=bernstein2(fr['coefficients'],*BOX)
                            assert all(v>=0 for row in bs for v in row)and any(v>0 for row in bs for v in row)
                        factor=mul(factor,power(fp,n))
                    assert mul(factor,q)==expected
                initial.append([bernstein2(p['coefficients'],*rect)for p in r['polynomials']])
            stats={'nodes':0,'leaves':0,'maxdepth':0}
            def walk(t,ws,depth,box):
                stats['nodes']+=1;stats['maxdepth']=max(stats['maxdepth'],depth)
                if type(t)is int:
                    assert 0<=t<len(ws)
                    if len(ws[t])==1:assert all(v>=0 for row in ws[t][0]for v in row)
                    else:assert witness_positive(ws[t],box)
                    stats['leaves']+=1;return
                if type(t)is dict:
                    check_combination_leaf(t,ws,box);stats['leaves']+=1;return
                assert type(t)is list and len(t)==3 and t[0]in(0,1)
                axis=t[0];ps=[[split2(b,axis)for b in w]for w in ws];l,r,bot,top=box
                left=(l,(l+r)/2,bot,top)if axis==0 else(l,r,bot,(bot+top)/2)
                right=((l+r)/2,r,bot,top)if axis==0 else(l,r,(bot+top)/2,top)
                walk(t[1],[[p[0]for p in w]for w in ps],depth+1,left)
                walk(t[2],[[p[1]for p in w]for w in ps],depth+1,right)
            walk(chart['tree'],initial,0,rect);assert stats==chart['statistics'];leaves+=stats['leaves']
    print('EXACTLY VERIFIED',len(seen),'blow-up patterns;',sum(len(c['charts'])for c in data['certificates']),
          'covering charts;',leaves,'Bernstein rectangles',flush=True)

_patterns=None
def reasons(w):
    global _patterns
    if _patterns is None:
        data=json.loads((ROOT/'data/blowup_certificates.json').read_text())
        _patterns={pat_key(r):i for i,r in enumerate(data['certificates'])}
    return ['blowup_pattern_'+str(_patterns[p])for p in patterns(w)if p in _patterns]

if __name__=='__main__':
    check_certificates();known=json.loads((ROOT/'data/certificates_n6.json').read_text())
    assert not any(reasons(w)for w in known)
    print('EXACTLY VERIFIED: no known witness excluded',flush=True)
