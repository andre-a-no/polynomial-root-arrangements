"""Exact finite Rolle-word enumeration and published necessary restrictions."""

def rolle_words(n):
    nodes=[(i,j) for i in range(n) for j in range(n-i)]
    idx={v:k for k,v in enumerate(nodes)}
    req=[0]*len(nodes)
    for i,j in nodes:
        if j: req[idx[i,j]] |= 1<<idx[i,j-1]
        if i:
            req[idx[i,j]] |= 1<<idx[i-1,j]
            req[idx[i-1,j+1]] |= 1<<idx[i,j]
    allmask=(1<<len(nodes))-1
    out=[]
    def rec(mask,s):
        if mask==allmask:
            out.append(s); return
        for k,(i,j) in enumerate(nodes):
            if not ((mask>>k)&1) and req[k]&mask==req[k]:
                rec(mask|(1<<k),s+str(i))
    rec(0,'')
    return out

def ks_ok(word,n):
    """Kostov-Shapiro 2002 Theorem 1 on every derivative and reflection."""
    pos=[[k for k,s in enumerate(word) if int(s)==i] for i in range(n)]
    for i in range(n-3):
        m=n-i
        if all(pos[i][j+1]<pos[i+2][j] for j in range(m-2)):
            if not pos[n-3][1]<pos[n-1][0]: return False
        if all(pos[i][j+1]>pos[i+2][j] for j in range(m-2)):
            if not pos[n-3][1]>pos[n-1][0]: return False
    return True

