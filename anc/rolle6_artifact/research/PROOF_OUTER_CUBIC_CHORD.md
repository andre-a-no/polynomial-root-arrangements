# Alternating values at the extreme cubic roots and the inner inflections

Established 2026-11-10. This rule excludes one further reflection pair.

Use the normalization and orientation from `PROOF_LEFT_MEAN_AVERAGE.md`:

\[
P=F+Bx+C,\quad F''=(x+u)x(x-1)(x-1-u-2r),\quad u,r>0,
\]
\[
H=u^2+u+2ur-2r-1.
\]

Let a<c<d be the roots of F'''. Assume H>0 and 0<r<1. Equivalently, P'' is negative at the left root of P^(4), and positive at the common mean mu=(1+r)/2.

**Lemma.** The four signs

\[
P(a)<0,\quad P(0)>0,\quad P(1)<0,\quad P(d)>0
\]

are impossible.

Define L(x)=F(x)-xF(1). Subtracting from P its secant through 0 and 1 shows that the proposed signs force L(a)<0<L(d), because a<0<1<d. Reduce L modulo F''', obtaining a quadratic A x²+B0 x+C0. Since a+d=3(1+r)/2-c and ad=c²-3(1+r)c/2-H/2, its values at a,d have the form

\[
L(a)=Va+W,\quad L(d)=Vd+W,
\]
\[
V=B_0+A[3(1+r)/2-c],\qquad
W=C_0-A[c^2-3(1+r)c/2-H/2].
\]

Thus a counterexample requires V>0, and since a>-u, also W-uV<0.

## Unbounded part: u>=2

Set D0=3+6r+2H. We have

\[
R'(1/2)=r/2>0,\quad R'(\mu)=-r(r+u)(r+u+1)<0,
\]

where R=F''. Therefore 1/2<c<mu. Since R'' is decreasing on (1/2,c),

\[
0=R'(c)\le r/2-D_0(c-1/2),\qquad
c\le1/2+r/(2D_0).
\]

V is affine in c. Exact coefficient identities show that both polynomials

\[
960\,\partial_c V,\qquad
-\frac{1920D_0}{r}V\left(u,r,\frac12+\frac{r}{2D_0}\right)
\]

have strictly positive coefficients after u=2+U. Their constant coefficients are positive. Hence V is increasing in c and is already negative at the upper bound on c. This contradicts V>0 throughout u>=2.

## Bounded part: 3/5<u<2

H>0 implies u>3/5. The equation R'(c)=0 gives

\[
D=u-2(u-1)c-3c^2,\quad
N=(2c-1)(u^2+u+2c-2c^2),\quad r=N/(2D).
\]

Because 1/2<c<1, N>0, so D>0. D decreases for c>=1/2, while D(u,3/5)=3/25-u/5<0. Consequently 1/2<c<3/5.

After substituting r=N/(2D), put

\[
\widehat V=30720D^5V,\quad \widehat W=30720D^5W.
\]

These are polynomials, and V-hat is divisible by 2c-1>0. A counterexample makes all six polynomials

\[
-D,\ -N,\ N-2D,\ -DH,\ -\widehat V/(2c-1),\ \widehat W-u\widehat V
\]

strictly negative. But the rectangle [3/5,2] x [1/2,3/5] has an exact **six-rectangle cover** on each leaf of which one of these polynomials has nonnegative Bernstein coefficients. Two leaves use N-2D, three use -V-hat/(2c-1), and one uses W-hat-uV-hat. The contradiction proves the lemma.

`outer_cubic_chord_obstruction.py` independently reconstructs the cubic and quadratic remainders, checks every cleared identity and the two tail polynomials, and computes each leaf's Bernstein bound both recursively and directly. All verification uses exact rational arithmetic and the standard library. The complete certificate is `data/outer_cubic_chord_certificate.json`.

The new pair is `012031402531024103210` and `012301420135204130210`. Neither it nor the other newly excluded pairs is a new realization. The combined continuation removes ten words from the preceding 38, leaving 28: 2706 <= a6 <= 2734. No exact count or publication novelty is claimed.
