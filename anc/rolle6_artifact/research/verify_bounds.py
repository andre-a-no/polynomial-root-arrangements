"""Independent obstruction checker: standard library and exact integers only.

Linear certificate: sum y_i E_i(g) = sum z_i B_i(g), z>=0, z!=0.
Quadratic certificate: a combination of E_i(g)*g_j and variance
identities equals a nonzero nonnegative combination of positive products.
All B_i(g) and all gaps are strictly positive in a realizable word.
"""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def forms(w,n):
    N=len(w);d=N-1;center=w.index(str(n-1))
    pos=[[k for k,c in enumerate(w) if int(c)==i] for i in range(n)]
    E=[[sum(k>j for k in pos[i])-(n-i)*int(center>j) for j in range(d)] for i in range(n-1)]
    B=[[int(i==j) for j in range(d)] for i in range(d)]
    def append_terms(terms):
        B.append([sum(a for k,a in terms if k>j) for j in range(d)])
    for i in range(n-2):
        m=n-i
        for k in range(1,m):
            append_terms([(pos[i+1][k-1],m-k+1),(pos[i][k-1],-(m-k)),(pos[i][k],-1)])
            append_terms([(pos[i][k-1],1),(pos[i][k],k),(pos[i+1][k-1],-(k+1))])
        for k in range(1,m-1):
            append_terms([(j,k) for j in pos[i][:k+1]]+[(j,-k-1) for j in pos[i+1][:k]])
            append_terms([(j,-k) for j in pos[i][-(k+1):]]+[(j,k+1) for j in pos[i+1][-k:]])
    return E,B

def check_linear(n):
    data=json.loads((ROOT/f'data/linear_rejected_n{n}.json').read_text())
    for w,v in data.items():
        E,B=forms(w,n);y=v[:len(E)];z=v[len(E):]
        assert len(v)==len(E)+len(B) and all(isinstance(x,int) for x in v)
        assert min(z)>=0 and sum(z)>0
        for j in range(len(w)-1):
            assert sum(a*e[j] for a,e in zip(y,E))==sum(a*b[j] for a,b in zip(z,B)),w
    print('EXACTLY VERIFIED',len(data),'linear obstruction certificates')
    return set(data)

def check_quadratic(n,pair_products=False):
    kind='quadratic_pairs' if pair_products else 'quadratic'
    data=json.loads((ROOT/f'data/{kind}_rejected_n{n}.json').read_text())
    for w,certificate in data.items():
        E,B=forms(w,n);d=len(w)-1
        pairs=[(i,j) for i in range(d) for j in range(i,d)]
        def product(a,b):
            return [a[i]*b[j]+(a[j]*b[i] if i!=j else 0) for i,j in pairs]
        I=B[:d];center=w.index(str(n-1))
        coords=[[int(j<k)-int(j<center) for j in range(d)] for k in range(len(w))]
        S=[]
        for i in range(n-1):
            roots=[product(coords[k],coords[k]) for k,c in enumerate(w) if int(c)==i]
            S.append([sum(r[j] for r in roots) for j in range(len(pairs))])
        numeq=len(E)*d+n-2
        total=[0]*len(pairs);positive_weight=0
        assert len({index for index,_ in certificate})==len(certificate)
        for index,weight in certificate:
            assert isinstance(index,int) and isinstance(weight,int) and index>=0
            if index<len(E)*d:
                i,j=divmod(index,d);row=product(E[i],I[j]);sign=1
            elif index<numeq:
                i=index-len(E)*d;m=n-i
                row=[2*a-m*(m-1)*b for a,b in zip(S[i],S[n-2])];sign=1
            else:
                assert weight>=0;positive_weight+=weight;sign=-1
                k=index-numeq
                if k<len(pairs):
                    i,j=pairs[k];row=product(I[i],I[j])
                elif k<len(pairs)+(len(B)-d)*d:
                    i,j=divmod(k-len(pairs),d)
                    assert i<len(B)-d
                    row=product(B[d+i],I[j])
                else:
                    assert pair_products
                    k-=len(pairs)+(len(B)-d)*d
                    bpairs=[(i,j) for i in range(d,len(B)) for j in range(i,len(B))]
                    assert k<len(bpairs)
                    i,j=bpairs[k];row=product(B[i],B[j])
            for j,x in enumerate(row):total[j]+=sign*weight*x
        assert positive_weight>0 and not any(total),w
    rejected=set(data)|{w[::-1] for w in data}
    print('EXACTLY VERIFIED',len(data),kind,'certificates, excluding',len(rejected),'words including reflections')
    return rejected

if __name__=='__main__':check_linear(6);check_quadratic(6)
