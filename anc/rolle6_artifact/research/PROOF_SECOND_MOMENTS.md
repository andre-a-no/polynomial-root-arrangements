# Coupling second moments of factors: two exact exclusions

8 November 2026. This stage starts with the 2442 unresolved words after the
metric stage. The objects counted are nondegenerate degree-six Rolle words,
not connected components of coefficient space. The moment identities and
Cauchy and triangle inequalities used here are elementary; no novelty claim
is made for them or for the resulting finite bounds.

## The exact moment identity

For a monic real-rooted polynomial f of degree m>=2, let mu be its root mean,
V(f)=sum_i (r_i-mu)^2, and d(f)>0 the distance between the two roots of
f^(m-2). Then

    4 V(f) = m(m-1) d(f)^2.

Indeed, translate mu to zero and write f=x^m+c2 x^(m-2)+.... Vieta gives
V(f)=-2c2. The monic normalization of f^(m-2) is
x^2+2c2/[m(m-1)], so its squared root gap is -8c2/[m(m-1)].
The identity is also valid for each derivative of each root-subset factor.
There are 72 occurrences of degree at least three in the shared model.

## Rational linear consequences

For an integer vector a with sum(a)=0, Cauchy's inequality gives

    2 sum_i a_i r_i <= sqrt(m(m-1) sum_i a_i^2) d(f).

Choose integer n with n^2>64^2 m(m-1)sum(a_i^2). Then
n*d(f)-128*sum(a_i*r_i)>0. In the implemented fixed family, take every
nonnegative integer vector of m-1 gaps summing to 4, form its cumulative
root coordinates x_i, and put a_i=m*x_i-sum(x). Divide by the gcd.
This gives 855 strict inequalities across the 72 occurrences.

An upper bound follows by writing the centered root vector as a positive
linear combination of centered step vectors. If g_k=r_(k+1)-r_k>0,

    d(f) <= sum_{k=1}^{m-1} [2 sqrt(k(m-k))/(m sqrt(m-1))] g_k.

This is the Euclidean triangle inequality. Round each coefficient strictly
upward with denominator 1024, using integer square roots. All 72 resulting
linear inequalities are strict. These 927 inequalities, on their own with
the previous constraints, gave no exclusions among the 2442 tested words.

This negative experiment does not prove that all possible second-moment
consequences are exhausted. An adaptive version on 66 sampled words also
found no exclusions; the four near-complete target families remained
unresolved in a deeper adaptive search. Feasibility in these relaxations
is not a polynomial realization certificate.

## Coupling a factor and a one-root deletion

Let S have n>=3 roots, T be S with root r deleted, mu_S its mean, and write
D=d(F_S), e=d(F_T), R=r-mu_S. Decomposing the centered sum of squares gives

    V(F_S) = V(F_T) + n/(n-1) R^2,

and therefore the exact identity

    D^2 = (n-2)/n e^2 + 4/(n-1)^2 R^2.

Put M=n(n-1)^2, alpha=(n-2)(n-1)^2, beta=4n. For any integers t,s,

    (alpha*t*e + beta*s*R)^2
       <= (alpha*t^2+beta*s^2) M D^2.

Thus, with N^2>64^2 M(alpha*t^2+beta*s^2),

    N*D - 64(alpha*t*e+beta*s*R) > 0.

The directions are (1,0), (1,s) for s=+-1,+-2,+-4,+-8, and (2,+-1).
There are 156 one-root deletions of subsets of size at least three. These
produce 1716 inequalities.

If r is the smallest or largest root of S, the sign of R is known strictly.
The triangle inequality also gives

    D <= sqrt((n-2)/n) e + 2/(n-1) |R|.

Round the square-root coefficient strictly upward with denominator 1024.
For the two extreme deletions of each of 42 subsets, this gives 84 more
strict linear inequalities. Total: 1800 coupling inequalities.

The complete model has 245 variables, 129 equalities, and 5892 strict
inequalities: 3165 from the preceding stage, 927 individual moment bounds,
and 1800 coupling bounds. All integer rounding conditions are checked by
the standard-library verifier.

## Exact result

The complete search rejected exactly these two candidates:

    012034102531042103210
    012301240135201430210

They are reflections. Each has an exact sparse integer Farkas certificate
in `data/moment_rejected_n6.json`, independently reconstructed and verified
by `research/moment_obstruction.py`. The nonnegative combination of strict
inequalities is nonempty and equals a combination of equalities, giving
0>0 for any putative realization. Floating-point LP is used only to find
the certificate; all final identities are checked coordinatewise with
Python integers.

Direct validation on 100 real-rooted polynomials checked 7200 moment
identities, 85500 lower and 7200 upper inequalities. Another independent
100-polynomial sample checked 15600 deletion identities and 180000 coupling
inequalities. These are implementation controls, not the proof.

Both excluded words were in the saved unresolved input and neither is
among the 2706 exact existence witnesses. This stage changes the upper
bound from 5148 to 5146. The following quartic-middle-root stage gives a
further improvement; its 82 certified exclusions are disjoint from these two.
