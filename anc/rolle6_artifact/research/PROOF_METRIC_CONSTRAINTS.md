# Metric restrictions on cubic chains and derivative spans

This stage strengthens the joint subfactor model of
`PROOF_SUBFACTOR_LATTICE.md`. It concerns words of all roots of a degree-six
strictly real-rooted polynomial and all its derivatives, with all 21 roots
distinct. It does not count connected components of coefficient space.

## 1. Exact cubic identity and rational cuts

Let a < b < c be the roots of a cubic, and u < v the roots of its derivative.
Put A=b-a, B=c-b, d=v-u. The quadratic formula gives

    9 d^2 = 4(A^2 + AB + B^2).

For positive integers t,s, put N=t^2+ts+s^2 and
L=(2t+s)A+(t+2s)B. Direct expansion gives

    4N(A^2+AB+B^2) - L^2 = 3(sA-tB)^2 >= 0.

Thus 3 sqrt(N) d >= L. Choose positive integers m,n with
n^2 > m^2 N. Because d>0,

    3n d - m(2t+s)A - m(t+2s)B > 0.

Every coefficient here is an integer. Strictness follows from the strict
rational rounding, even when the underlying non-strict inequality is an
equality. The implemented choices are (t,s,m)=(1,1,56),
(1,k,64),(k,1,64) for k=2,4,8,16, and (3,5,64),(5,3,64).
For each choice n=floor(m sqrt(N))+1 is computed by integer square root;
coefficients are divided by their common divisor.

For example t=s=1 gives

    97(v-u) > 56(c-a).

All 11 inequalities are applied to every cubic in every derivative chain of
every root-subset factor of p. There are 42 such cubics: p''' (one), the second
derivatives of its six quintic factors, the first derivatives of its fifteen
quartic factors, and its twenty cubic factors. This adds 462 necessary strict
linear inequalities to the 3093 inequalities and 129 equalities of the
previous model, in the same 245 variables.

## 2. A classical span bound, also tested

For a real-rooted polynomial f of degree m>=3, let span(f) be its largest
root minus its smallest root, with multiplicities allowed. Then

    span(f') >= sqrt((m-2)/m) span(f).

This is a special case of the classical theorem of A. Meir and A. Sharma,
“On Zeros of Derivatives of Polynomials”, Canadian Mathematical Bulletin
11(3) (1968), 443–445, DOI https://doi.org/10.4153/CMB-1968-051-x.
The theorem treats all derivative orders. The publisher's original PDF and
its first-page theorem were checked. Neither this span inequality nor the
cubic quadratic formula is claimed as new.

Here is a self-contained verification of the first-derivative bound.
Normalize the extreme roots to 0 and 1 and vary the m-2 remaining labeled
roots over [0,1]^(m-2). The span of the derivative is continuous and attains
a minimum. If an interior root reaches 0, the smallest derivative root is 0;
the largest derivative root v satisfies v>=(m-1)/m. Indeed if 1 is simple,
the logarithmic-derivative equation gives

    1/(1-v) = sum_{i=1}^{m-1} 1/(v-x_i) >= (m-1)/v.

If 1 is repeated, v=1. The boundary with a root at 1 is treated by reflection.
Thus every boundary point has derivative span at least (m-1)/m, larger than
sqrt((m-2)/m). A smaller value is realized by putting every interior root at
1/2, so the minimum is attained in the open cube.

Write u,v for the extreme critical roots and S_z=sum_j (z-x_j)^(-2).
Both u and v are simple, even when interior roots coincide. Implicit
differentiation gives du/dx_i=(u-x_i)^(-2)/S_u and
dv/dx_i=(v-x_i)^(-2)/S_v. At a minimum their difference is zero for each
interior root. Hence (v-x_i)/(x_i-u)=sqrt(S_u/S_v), independent of i.
This strictly decreasing function of x_i forces all interior roots to be
equal to t. For f=x(x-1)(x-t)^(m-2), the extreme critical roots are the roots
of m x^2-(m-1+2t)x+t; their squared gap is

    [m(m-2)+4(t-1/2)^2]/m^2.

Its minimum is (m-2)/m, as asserted. For m>3, equality requires repeated
interior roots. The rational cuts below are strict regardless of equality.

For m=4,5,6, respectively, use the rational lower approximations
70/99, 189/244, 218/267. Their validity is checked using integer squares.
If q denotes the denominator and r the numerator, impose

    q span(f') - r span(f) > 0.

These add 30 rows (22 quartic, 7 quintic, 1 sextic), for a total of 3585
strict inequalities. The cubic span row is already among the cubic cuts.

The simpler model keeps only one span inequality per derivative link:
the first of the 11 cubic cuts for each of 42 cubics, and the 30 higher-degree
span cuts. It has 72 new rows, hence 3165 strict inequalities in total.
The expanded model with all cubic cuts has 3585. Both are certified below.

## 3. Exact infeasibility certificates

The numerical search solves the homogeneous relaxation E z=0, B z>0 by
rescaling to B z>=1. This is equivalent because there are finitely many
strict inequalities. Numerical infeasibility alone is not accepted.

Each stored certificate consists of integer equality multipliers y_i
(of either sign) and positive integer inequality multipliers c_j, with
nonempty inequality support, satisfying exactly

    sum_i y_i E_i = sum_j c_j B_j.

Every real realization would make the left-hand side zero and the right-hand
side strictly positive. The identity is checked coordinate by coordinate
with Python integers in `research/metric_obstruction.py`. The exact verifier
reconstructs all matrices using only the standard library and is independent
of SciPy and SymPy used in the search and rational reconstruction.

The proof requires the theoretical necessity of the constraints above and
in `PROOF_SUBFACTOR_LATTICE.md`; a Farkas identity alone does not validate
an incorrectly chosen constraint.

## 4. Scope of the experiment

The input is the saved list of 2448 unresolved words after the subfactor-sign
stage. Files `research/work/cubic_metric_output.txt` and
`research/work/span_metric_output.txt` record the two searches. Only exact
certificates in `data/metric_rejected_n6.json` enter the bound.

The cubic identity itself is nonlinear. The eleven rational inequalities
are only necessary consequences of it; feasibility of the relaxation is
not a realization proof. Increasing the number of directions is not known
to resolve every infeasible word. No claim of publication novelty or of
logical independence from all previous approaches is made. The set audit
checks that these stored exclusions were not previously counted.

## 5. Certified result and overlap

The cubic-only experiment rejects four words. The expanded experiment
rejects six. All six have also been certified using only the 72 span cuts;
thus the additional cubic directions make no additional exclusions in this
2448-word input. The exact data files `metric_rejected_n6.json` (span model)
and `metric_full_rejected_n6.json` (expanded model) contain the same six words.
`check_certificates()` verifies both sets and counts their union only once:

    010234123012051403210
    012034102531402103210
    012034123012051403210
    012301204135201430210
    012304150210321430210
    012304150210321432010

The set is invariant under reflection. It is disjoint from all 2706 exact
existence witnesses and from the earlier exclusions. It changes the proven
interval from [2706,5154] to [2706,5148], leaving 2442 unresolved words.
No new derivative-word/mean-sign family is completely classified; there
remain 16 fully classified families with 254 realized words, already included
in 2706. The four previously discussed families with 17 realized words and
one unresolved word each are still unresolved.

The full expanded search reported zero solver-unknown outcomes. In the
six-word overlap probe, the span-only model rejected all six; since it is a
subset of the expanded model, no extra span-only exclusion is expected among
expanded-feasible words. This is a statement about the numerical search,
not a mathematical feasibility certificate for the remaining 2442 words.
The certified upper bound needs only the six exact contradictions.
