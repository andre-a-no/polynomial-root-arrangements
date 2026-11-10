"""Exact sign deductions from every product splitting of every subproduct.
Coefficients in Leibniz's formula are positive, so only signs are needed here.
"""
import itertools
FULL=tuple(range(6))
SUBSETS=[s for n in range(1,7)for s in itertools.combinations(FULL,n)]
NODES=[(S,j)for S in SUBSETS for j in range(len(S)+1)]
INDEX={v:i for i,v in enumerate(NODES)}
PRODUCTS=[]
for U in SUBSETS:
    n=len(U)
    if n<2:continue
    for m in range(1,1<<n):
        S=tuple(U[i]for i in range(n)if m>>i&1)
        T=tuple(x for x in U if x not in S)
        if not T or S>T:continue
        for j in range(1,n):
            terms=[(INDEX[S,a],INDEX[T,j-a])for a in range(j+1)if a<=len(S)and j-a<=len(T)]
            PRODUCTS.append((INDEX[U,j],terms))


def propagate(w):
    pos=[[i for i,c in enumerate(w)if c==str(j)]for j in range(6)]
    signs=[[None]*21 for v in NODES];bounds={}
    for v,(S,j)in enumerate(NODES):
        n=len(S)
        if j==n:signs[v]=[1]*21
        elif j==0:
            points=[pos[0][i]for i in S]
            signs[v]=[0 if x in points else (-1)**sum(z>x for z in points)for x in range(21)]
            bounds[v]=[[z,z]for z in points]
        elif S==FULL:
            signs[v]=[0 if x in pos[j]else (-1)**sum(z>x for z in pos[j])for x in range(21)]
            bounds[v]=[[z,z]for z in pos[j]]
        else:bounds[v]=[[pos[0][S[i]],pos[0][S[i+j]]]for i in range(n-j)]
    # Strict ordering constraints between roots. Original/subproduct roots are fixed.
    relations=[]
    for S in SUBSETS:
        n=len(S)
        for j in range(n):
            v=INDEX[S,j]
            for i in range(n-j-1):relations.append(((v,i),(v,i+1)))
        for j in range(n-1):
            a,b=INDEX[S,j],INDEX[S,j+1]
            for i in range(n-j-1):relations.extend([((a,i),(b,i)),((b,i),(a,i+1))])
        if n<3:continue
        for r in S:
            T=tuple(x for x in S if x!=r)
            for j in range(1,n-1):
                a,b=INDEX[S,j],INDEX[T,j]
                for i in range(n-j-1):relations.extend([((a,i),(b,i)),((b,i),(a,i+1))])
    for S in SUBSETS:
        for x in S:
            if x==5 or x+1 in S:continue
            T=tuple(sorted((set(S)-{x})|{x+1}))
            for j in range(1,len(S)):
                for i in range(len(S)-j):relations.append(((INDEX[S,j],i),(INDEX[T,j],i)))
    changes=0;rounds=0
    class Conflict(Exception):pass
    def set_sign(v,x,value):
        nonlocal changes
        old=signs[v][x]
        if old is None:signs[v][x]=value;changes+=1
        elif old!=value:raise Conflict(('sign',NODES[v],x,old,value))
    def bound(v,i,lo=None,hi=None):
        nonlocal changes
        b=bounds[v][i];fixed=b[0]==b[1]
        if fixed:
            if lo is not None and b[0]<=lo:raise Conflict(('fixed_lower',NODES[v],i,lo,b))
            if hi is not None and b[0]>=hi:raise Conflict(('fixed_upper',NODES[v],i,hi,b))
            return
        nl=b[0]if lo is None else max(b[0],lo)
        nh=b[1]if hi is None else min(b[1],hi)
        if nl>=nh:raise Conflict(('empty_interval',NODES[v],i,nl,nh))
        if b!=[nl,nh]:b[:]=[nl,nh];changes+=1
    try:
        while True:
            before=changes;rounds+=1
            for (va,ia),(vb,ib)in relations:
                a,b=bounds[va][ia],bounds[vb][ib]
                if a[0]==a[1]:
                    if b[0]==b[1]:
                        if not a[0]<b[0]:raise Conflict(('fixed_order',va,ia,vb,ib))
                    else:bound(vb,ib,lo=a[0])
                elif b[0]==b[1]:bound(va,ia,hi=b[0])
                else:
                    bound(va,ia,hi=b[1]);bound(vb,ib,lo=a[0])
            for v,row in bounds.items():
                n=len(row)
                for x in range(21):
                    if any(lo==hi==x for lo,hi in row):set_sign(v,x,0);continue
                    mn=sum(hi<x if lo==hi else hi<=x for lo,hi in row)
                    mx=sum(lo<x for lo,hi in row)
                    s=signs[v][x]
                    if mn==mx:set_sign(v,x,(-1)**(n-mn));continue
                    if s is None or s==0:continue
                    possible=[k for k in range(mn,mx+1)if (-1)**(n-k)==s]
                    if not possible:raise Conflict(('parity',NODES[v],x,mn,mx,s))
                    for i in range(min(possible)):bound(v,i,hi=x)
                    for i in range(max(possible),n):bound(v,i,lo=x)
            for parent,terms in PRODUCTS:
                for x in range(21):
                    unknown=[];nonzero=set()
                    for a,b in terms:
                        sa,sb=signs[a][x],signs[b][x]
                        if sa==0 or sb==0:continue
                        if sa is None or sb is None:unknown.append((a,b))
                        else:nonzero.add(sa*sb)
                    if len(nonzero)>1:continue
                    t=next(iter(nonzero))if nonzero else 0
                    if not unknown:set_sign(parent,x,t);continue
                    if len(unknown)!=1:continue
                    s=signs[parent][x]
                    if s is None:continue
                    if t==0:u=s
                    elif s==0:u=-t
                    elif s!=t:u=s
                    else:continue
                    a,b=unknown[0];sa,sb=signs[a][x],signs[b][x]
                    if sa is not None:set_sign(b,x,u*sa)
                    elif sb is not None:set_sign(a,x,u*sb)
            if changes==before:break
            assert rounds<100
    except Conflict as e:return {'conflict':e.args[0],'rounds':rounds}
    return {'bounds':{NODES[v]:row for v,row in bounds.items()},'signs':{NODES[v]:s for v,s in enumerate(signs)},'rounds':rounds,'deductions':changes}

if __name__=='__main__':
    import sys,json,time
    from pathlib import Path
    root=Path(__file__).resolve().parents[1];words=json.loads((root/'data/certificates_n6.json').read_text());start=time.time()
    for w in list(words)[::100]:
        r=propagate(w);assert 'conflict'not in r,(w,r)
    print('Controls passed',len(list(words)[::100]),'seconds',time.time()-start)
