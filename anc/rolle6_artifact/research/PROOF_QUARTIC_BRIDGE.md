# A quartic's middle critical root: using an earlier lemma throughout the factor model

8 November 2026. This stage uses the same nondegenerate degree-six Rolle-word
problem and the same exact existence witnesses as the previous stages.
Neither a new classical theorem nor publication priority is claimed.

## Elementary lemma

Let f have four simple real roots r1<r2<r3<r4. Let u<v<w be the roots of f',
mu=(r1+r2+r3+r4)/4 and kappa=(r2+r3)/2.
If kappa=mu, then v=kappa=mu. Otherwise v lies strictly between mu and kappa.

Proof. Scale the leading coefficient to one and translate mu to zero. Write
the translated roots as a<b<c<d, a+b+c+d=0. Direct expansion gives

    f'(0)=(a+b)(a+c)(b+c).

Both a+b and a+c are negative: their complementary sums are larger and
opposite. Thus sign(f'(mu))=sign(kappa-mu). Also, differentiating the product
of the outer-root and inner-root quadratic factors gives

    f'(kappa)=(2*kappa-r1-r4)(kappa-r2)(kappa-r3),

whose sign is -sign(kappa-mu). The mean mu is also the mean of u,v,w, so
u<mu<w. Interlacing gives u<r2<kappa<r3<w. Throughout (u,w), the sign of
f'(x) is the sign of v-x. Hence, if kappa>mu, then mu<v<kappa; if kappa<mu,
then kappa<v<mu. If they are equal, f'(mu)=0 and the zero must be v.

The two polynomial identities are checked exactly by the standard-library
verifier. A separate implementation control checked 10000 random quartics.

## Linear consequences when the side of the mean is known

Since 4(kappa-mu)=(r2-r1)-(r4-r3), the following hold strictly:

    v>mu => r2-r1 > r4-r3 and r2+r3 > 2v;
    v<mu => r2-r1 < r4-r3 and r2+r3 < 2v.

For the quartic p'', the original word specifies the order of its middle
critical root (the middle root of p''') and mu (the root of p^(5)).
They are distinct by the definition of the words being counted.

The shared model also contains the six quartics obtained by differentiating
quintic factors and the fifteen quartic factors themselves. For each of
these 21 auxiliary quartics we apply the lemma only when the previously
proved interval bounds put v strictly to one side of its own mean. If this
order is unknown, no strict inequality is added. In particular, coincidence
of auxiliary roots is not silently ruled out.

The interval rule is exact: if v is in (P_a,P_b) and mu is in (P_c,P_d),
then a>=d implies v>mu, and b<=c implies v<mu. The intervals come from the
existing subfactor sign proof and their linear inequalities remain in the
model. There is no inference from approximate coordinates here.

## Relation to the earlier work

The same location property was already implicit in the normalized
middle-interval proof in `PROOF_MIDDLE_AND_CONVEXITY.md`: that proof showed
that xi-nu has the opposite sign to xi-1/2, with equality simultaneously.
The earlier exclusion rule used signs of values of p and p'. Here we use
its quantitative root-location consequence as linear inequalities and
apply it to the shared family of factors as well. This is additional use
of an established lemma, not a logically independent theorem.

## Certificates and accounting

On the 2442-word unresolved input after the metric stage, this model found
82 exclusions and no solver-unknown outcomes. Every exclusion is stored as
an integer Farkas certificate in `data/quartic_bridge_rejected_n6.json`.
`research/quartic_bridge_obstruction.py` reconstructs the matrices without
NumPy or SciPy and checks the identities exactly. Its inputs include all
previous shared-factor constraints and the 72 span inequalities.

The two exclusions from the coupled-second-moment experiment are not among
these 82. The union therefore contains 84 new exclusions. No word previously
excluded and no exact existence witness is counted again. The resulting
bounds are

    2706 <= a6 <= 5064,

with 2358 unresolved words. The full check in `research/check_all.py`
recomputes the bounds, verifies all existence and nonexistence certificates,
and audits overlaps. The full list of closed families is recomputed from
the remaining words, rather than inferred from sample counts.

The standalone skewness-only sample used only the first of the two
consequences and found one exclusion among 66 sampled words. It is a
weaker part of the quartic lemma and is not counted as another method or
additional exclusion.

## Additional implementation controls

All 2706 words with exact existence witnesses passed the numerical
quartic-constraint model. Inspection of the stored certificate supports
shows that 79 of the 82 certificates use added rows for auxiliary quartics;
three use added quartic rows only for p''. This describes the certificates
found, not a theorem that every such use is indispensable.
The final full exact verification completed successfully, including the
set-based overlap audit and the unchanged 16 complete families.
