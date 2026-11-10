# A left cubic root and the common mean bound an inflection secant

Established 2026-11-10. This is an exact computer-assisted proof, not a numerical search exclusion. It removes four of the preceding 38 unresolved words.

Let the leading coefficient of P be positive. Write Q=P', R=P'', and let z1<z2<z3<z4 be the roots of R. Reflect if necessary so that the middle root of R' lies to the left of the common mean mu. All roots in the Rolle word are distinct, so this choice is strict. Let a be the leftmost root of R' and v- the left root of R''.

**Lemma.** If R(v-)<0 and R(mu)<0, then

\[
\frac{P(z_3)-P(z_1)}{z_3-z_1}
\leq\max\{Q(a),Q(\mu)\}.
\]

Consequently P(z1)<0<P(z3), Q(a)<0 and Q(mu)<0 cannot hold simultaneously with these hypotheses.

## Normalization and the physical domain

A positive affine change and positive rescaling give

\[
z_1=-u,\quad z_2=0,\quad z_3=1,\quad z_4=T=1+u+2r,
\qquad u,r>0.
\]

The orientation r>0 follows from the sign of the centered cubic moment, r(r+u)(r+u+1); this is also the orientation used in the earlier inflection proofs. Put

\[
H=u^2+u+2ur-2r-1,\qquad \mu=(1+r)/2,
\]
\[
F(x)=x^6/30-(1+r)x^5/10-Hx^4/12+uTx^3/6,
\quad P(x)=F(x)+Bx+C.
\]

Thus F''=R. R(v-)<0 is equivalent to H>0: the negative left root of R'' then lies in (-u,0), where R is negative. More explicitly R''=12x^2-12(1+r)x-2H, and R''(-u)=12u^2+12(1+r)u-2H>0, so its negative root is greater than -u. Conversely v- is less than the middle root of R', which lies in (0,1); when v- is nonnegative R(v-) is positive. Also R(mu)<0 is equivalent to r>1, since mu is between 0 and T. H>0 implies u>3/5.

Write a=-ut. The logarithmic derivative

\[
\frac{R'(-ut)}{R(-ut)}
=\frac1{u(1-t)}-\frac1{ut}-\frac1{ut+1}-\frac1{ut+T}
\]

is negative at t=1/2 and positive at t=3/4. There is exactly one root a in (-u,0), so 1/2<t<3/4. At that root,

\[
D=ut(3t-2)+2t-1
=\frac{ut(1-t)(ut+1)}{ut+T}>0,
\]
\[
N=(2ut+1)[u(1-2t^2)-2t+1],\qquad r=N/(2D).
\]

Set A0=(F(1)-F(-u))/(1+u). The constants B,C cancel in the desired inequality.

## Polynomial certificate

The certificate contains five polynomials, with definitions that do not depend on a symbolic algebra system:

\[
-D,\quad 2D-N,\quad -DH/u,\quad L,\quad M,
\]

where, after r=N/(2D),

\[
60D(F'(-ut)-A_0)=(2ut-u+1)L,
\]
\[
1920D^4(F'(\mu)-A_0)=(r+u)M.
\]

The factors D, u, 2ut-u+1 and r+u are positive. Therefore at a counterexample all five displayed polynomials would be strictly negative: the first three by the physical domain, the last two by failure of the lemma.

For reference,

\[
120(F'(\mu)-A_0)/(r+u)
=-3r^4-7r^3u-5r^3+2r^2u^2+5r^2-2ru^3+25ru+15r
-6u^4-12u^3+2u^2+8u+2.
\]

Compactify the unbounded parameter by

\[
u=\frac{3+2v}{5(1-v)},\qquad 0<v<1.
\]

For each polynomial of u-degree n, multiply its substitution by (1-v)^n. This preserves the sign throughout the physical domain and produces a polynomial on [0,1] x [1/2,3/4]. The closed rectangle, including all limiting faces, is covered by **27 rational rectangles**. On each, one of the five transformed polynomials has nonnegative Bernstein coefficients. The leaf counts by polynomial are **2,3,8,13,1**. This precludes their simultaneous strict negativity.

`left_mean_average_obstruction.py` reconstructs the cleared identities using rational polynomial arithmetic, verifies compactification coefficient by coefficient, and checks each leaf both by recursive Bernstein subdivision and an independent direct Bernstein formula. It uses only the Python standard library. The complete coefficients and split tree are in `data/left_mean_average_certificate.json`.

## Application and overlap

The rule matches 130 of the 13404 initial candidates. Of these, 126 were already excluded; four are new:

- `012010342130210543210`
- `012013042130210543210`
- `012345012031240310210`
- `012345012031243010210`

They form two reflection pairs. None of the 2706 exact realizations satisfies the forbidden conditions. After this rule alone the bounds are 2706 <= a6 <= 2740 and 34 words remain. Subsequent rules are reported separately.
