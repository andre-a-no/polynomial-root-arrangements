"""Exact obstruction to a nonnegative stationary value in (0,1).

Six polynomial identities and 65 rational Bernstein rectangles. The
verifier uses only Python's standard library; SymPy is not imported.
"""
from fractions import Fraction as F
from pathlib import Path
import json
from two_parameter_obstruction import add, scale, mul, power, dense_to_sparse, bernstein2, split2

ROOT = Path(__file__).resolve().parents[1]
NAMES = ['Tc-Vl', 'Tc-V0', 'Tl-V0', 'Tl-Ul', 'Tc-Ul', 'Tc-Ua']


def generated_comparisons():
    one = {(0, 0): 1}
    K, U = {(1, 0): 1}, {(0, 1): 1}
    D = add(power(K, 2), scale(K, -1), one)
    AN = add(scale(one, 2), scale(K, -2), scale(power(K, 2), -1))
    CN = add(scale(one, -1), scale(K, -2), scale(power(K, 2), 2))
    DN = scale(add(AN, CN), -1)
    BN = add(scale(mul(CN, power(D, 2)), 3), scale(power(CN, 3), -1))
    assert not add(power(AN, 2), mul(AN, CN), power(CN, 2), scale(power(D, 2), -3))
    assert not add(mul(mul(AN, CN), DN), BN)
    X = mul(U, D)
    L = scale(D, -1)

    def sums(Y, n, weighted=False):
        return add(*(scale(mul(power(Y, n-i), power(X, i)), i+1 if weighted else 1) for i in range(n+1)))

    # These are 60 D^4 times the quantities in the proof.
    def T(Y):
        return add(scale(sums(Y, 4), -3), scale(mul(power(D, 2), sums(Y, 2)), 30),
                   scale(mul(BN, add(Y, X)), -30))

    def Upper(Y):
        return add(scale(sums(Y, 4, True), -1), scale(mul(power(D, 2), sums(Y, 2, True)), 15),
                   scale(mul(BN, add(Y, scale(X, 2))), -20))

    def V(Y):
        return add(scale(power(Y, 4), -15), scale(mul(power(D, 2), power(Y, 2)), 90),
                   scale(mul(BN, Y), -60))

    # Independently verify the divided-difference formulas against P0.
    def P120(Y):
        return add(power(Y, 6), scale(mul(power(D, 2), power(Y, 4)), -15),
                   scale(mul(BN, power(Y, 3)), 20))

    def Q60(Y):
        return add(scale(power(Y, 5), 3), scale(mul(power(D, 2), power(Y, 3)), -30),
                   scale(mul(BN, power(Y, 2)), 30))

    for Y in (AN, L, CN, {}):
        delta = add(Y, scale(X, -1))
        assert mul(T(Y), delta) == add(Q60(X), scale(Q60(Y), -1))
        assert mul(Upper(Y), power(delta, 2)) == add(
            scale(P120(Y), -1), P120(X), scale(mul(Q60(X), delta), 2))

    comparisons = [add(T(CN), scale(V(L), -1)), add(T(CN), scale(V({}), -1)),
                   add(T(L), scale(V({}), -1)), add(T(L), scale(Upper(L), -1)),
                   add(T(CN), scale(Upper(L), -1)), add(T(CN), scale(Upper(AN), -1))]
    factors = [scale(one, 3), scale(one, 3), scale(D, 3), mul(D, add(one, U)), one, one]
    return comparisons, factors


def reasons(word):
    w = word
    if [i for i, c in enumerate(w) if c == '3'][1] > w.index('5'):
        w = w[::-1]
    pos = [[i for i, c in enumerate(w) if c == str(j)] for j in range(6)]
    a, c, d = pos[3]
    left, right = pos[4]
    mu = pos[5][0]
    assert a < left < c < mu < right < d

    def sign(j, x):
        return (-1) ** sum(i > x for i in pos[j])

    required = ((0, a), (0, left), (1, left), (1, c), (2, left), (2, mu))
    if any(sign(j, x) != -1 for j, x in required):
        return []
    if not any(mu < x < right and sign(0, x) > 0 for x in pos[1]):
        return []
    return ['right_positive_stationary_value']


def check_certificates():
    data = json.loads((ROOT / 'data/right_maximum_certificates.json').read_text())
    records = data['polynomials']
    assert [r['id'] for r in records] == NAMES
    comparisons, factors = generated_comparisons()
    for record, comparison, factor in zip(records, comparisons, factors):
        assert all(type(v) is int for row in record['coefficients'] for v in row)
        assert mul(factor, dense_to_sparse(record['coefficients'])) == comparison, record['id']
    assert data['rectangle'] == [[1, 1], [11, 8], [0, 1], [1, 1]]
    rect = [F(*v) for v in data['rectangle']]
    assert 3 < (2 * rect[1] - 1) ** 2  # (1+sqrt(3))/2 < 11/8
    initial = [bernstein2(r['coefficients'], *rect) for r in records]
    counts = dict.fromkeys(NAMES, 0)
    nodes, depthmax = 0, 0

    def walk(tree, bs, depth):
        nonlocal nodes, depthmax
        nodes += 1
        depthmax = max(depthmax, depth)
        if type(tree) is int:
            assert 0 <= tree < len(records)
            assert all(v >= 0 for row in bs[tree] for v in row)
            counts[NAMES[tree]] += 1
            return
        assert isinstance(tree, list) and len(tree) == 3
        axis, left, right = tree
        assert type(axis) is int and axis in (0, 1)
        pairs = [split2(b, axis) for b in bs]
        walk(left, [p[0] for p in pairs], depth + 1)
        walk(right, [p[1] for p in pairs], depth + 1)

    walk(data['tree'], initial, 0)
    assert data['statistics'] == {'nodes': nodes, 'leaves': sum(counts.values()),
                                 'maxdepth': depthmax, 'by_polynomial': counts}
    print('EXACTLY VERIFIED 6 right-maximum identities and', sum(counts.values()),
          'Bernstein rectangles covering the full domain', flush=True)


if __name__ == '__main__':
    check_certificates()
    known = json.loads((ROOT / 'data/certificates_n6.json').read_text())
    assert not any(reasons(w) for w in known)
    print('EXACTLY VERIFIED: no known witness is excluded', flush=True)
