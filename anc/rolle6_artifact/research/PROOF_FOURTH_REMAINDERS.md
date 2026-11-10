# Two obstructions from the remainder modulo the fourth derivative

Established 2026-11-10. These are exact sign arguments with rational polynomial certificates. Together they exclude four additional words, beyond the four excluded by the left-cubic/mean secant lemma.

Use the same orientation and normalization

\[
R=P''=(x+u)x(x-1)(x-T),\quad T=1+u+2r,\quad u,r>0,
\]
\[
H=u^2+u+2ur-2r-1,\quad
F=x^6/30-(1+r)x^5/10-Hx^4/12+uTx^3/6,\quad P=F+Bx+C.
\]

Let v-<v+ be the roots of P^(4)=R''. Assume R(v-)<0, equivalently H>0. Then u>3/5 and

\[
v_\pm=\mu\pm\delta,\quad
\mu=(1+r)/2,\quad \delta^2=\mu^2+H/6,
\quad v_-<0<1<v_+.
\]

The last inequality follows from R''(1)=-12r-2H<0. At either root of R'',

\[
x^2=(1+r)x+H/6.
\]

All remainders below are obtained by this quadratic identity; the verifier reconstructs them rather than trusting symbolic output.

## 1. The derivative average

**Lemma.** Under H>0,

\[
P(1)-P(0)>\min\{P'(v_-),P'(v_+)\}.
\]

The interval has length 1; equivalently its mean derivative exceeds at least one of these two values. Hence P(0)>0>P(1) together with P'(v-)>0 and P'(v+)>0 is impossible.

Write F'(v)=A v+B0 at the two roots of R''. Then

\[
F(1)-\frac{F'(v_-)+F'(v_+)}2=\frac{rS}{60},
\]
\[
S=9r^4+30r^3u+15r^3+35r^2u^2+35r^2u+5r^2
+20ru^3+30ru^2-5r+5u^4+10u^3-5u-1.
\]

For u>=7/10 all coefficients of S(7/10+U,r) are positive, including the constant 261/2000. Thus the mean of the two endpoint values is strictly smaller than F(1), proving the lemma in that range.

For 3/5<u<7/10, H>0 implies

\[
0<r<\frac{u^2+u-1}{2(1-u)}<19/60<1/3.
\]

On the entire rectangle [3/5,7/10] x [0,1/3], exact Bernstein coefficients give

\[
20A\ge882373/270000>0,\qquad
120(F(1)-B_0)\ge7347/625>0.
\]

Since v-<0, F'(v-)=A v-+B0<B0<F(1). Adding the same B to the derivative and the secant proves the claim for P.

This rule matches 66 initial candidates, of which 64 were already excluded. The new pair is

- `012301402531021403210`
- `012304120135204103210`.

## 2. Alternating values at two fourth-derivative roots and two inflections

**Lemma.** If H>0 and 0<r<1, the signs

\[
P(v_-)<0,\quad P(0)>0,\quad P(1)<0,\quad P(v_+)>0
\]

are impossible. The condition r<1 is exactly R(mu)>0 in the chosen orientation.

Subtract from P its line through 0 and 1. Affine terms disappear. The four proposed signs imply

\[
F(v_-)-v_-F(1)<0,\qquad F(v_+)-v_+F(1)>0.
\]

Indeed the corresponding expression for P is P(v)-(1-v)P(0)-vP(1); each term has the indicated strict sign at v-<0 and v+>1.

Reduce F modulo R'' and write F(v)=A1 v+B1. Put d=A1-F(1). Opposite signs as above force d>0. Direct identities give

\[
d=-rS_1/180,
\]
\[
S_1=12r^4+40r^3u+20r^3+40r^2u^2+40r^2u
+20ru^3+30ru^2-30ru-20r+5u^4+10u^3-15u^2-20u-3,
\]
\[
J=3240(1+r)B_1/H-540d.
\]

All coefficients of S1(2+U,r) are positive. Thus u>=2 gives d<0, a contradiction. For 3/5<u<2 and 0<r<1, on the full enclosing rectangle exact Bernstein coefficients give

\[
J\ge32843/625>0.
\]

Finally

\[
-v_-=\delta-\mu=\frac{H}{6(\delta+\mu)}
<\frac{H}{6(1+r)}.
\]

If d>0, it follows that

\[
F(v_-)-v_-F(1)=B_1+dv_-
>B_1-\frac{dH}{6(1+r)}
=\frac{HJ}{3240(1+r)}>0,
\]

contradicting the first necessary inequality.

This rule matches 30 initial candidates, of which 28 were already excluded. Its new pair is

- `012304102531021043210`
- `012340120135201403210`.

## Verification

`fourth_remainder_obstruction.py` reconstructs the quadratic remainders of F and F', checks the displayed identities, checks the two positive-coefficient tails, and independently computes all three Bernstein lower bounds in two ways on the two stated rectangles. It uses only exact rational arithmetic from the standard library.

The coefficients are stored in `data/fourth_remainder_certificate.json`. Neither rule matches any of the 2706 exact realizations. The new pairs and the four words of the preceding lemma are disjoint. Combined: eight new exclusions, 30 unresolved words, 2706 <= a6 <= 2736. This remains an interval, not an exact enumeration theorem.
