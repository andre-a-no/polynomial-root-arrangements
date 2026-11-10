"""Exact four-rectangle proof of a first-minimum obstruction.

No numerical roots, optimizers, or symbolic algebra packages are required.
"""
import json
from pathlib import Path
from fractions import Fraction as F
from two_parameter_obstruction import add, scale, mul, power, bernstein2
from stationary_obstruction import dense
from marker_hermite_obstruction import decode
from marker_hermite6_obstruction import univariate_division
from inflection_obstruction import patterns

ROOT = Path(__file__).resolve().parents[1]
MASK = [[1, -1], [3, 1], [6, -1], [7, 1], [11, -1]]


def polynomials():
    one = {(0, 0): 1}
    k, z = {(1, 0): 1}, {(0, 1): 1}
    km = add(k, scale(one, -1))
    D = add(power(k, 2), scale(k, -1), one)
    zp = add(z, one)
    U = add(scale(mul(power(k, 2), zp), 30),
            scale(mul(mul(D, zp), add(power(z, 2), one)), -5),
            scale(mul(km, add(*(power(z, j) for j in range(5)))), -3))
    zm = add(z, scale(k, -1))
    L = add(scale(mul(power(k, 2), zm), 30),
            scale(mul(mul(D, zm), add(power(z, 2), power(k, 2))), -5),
            scale(mul(km, add(*(scale(mul(power(z, 4-j), power(k, j)), (-1)**j)
                                for j in range(5)))), -3))
    M = add(scale(power(k, 2), 40), scale(mul(D, power(z, 2)), -8),
            scale(mul(km, power(z, 3)), -5))
    Q, rem = univariate_division(scale(add(U, scale(L, -1)), -1), add(k, one), 0)
    assert not rem
    numerator = add(scale(add(mul(z, M), scale(U, -1)), 10),
                    scale(mul(add(scale(z, 3), scale(one, -5)), add(U, scale(L, -1))), -1))
    N, rem = univariate_division(numerator, km, 0)
    assert not rem
    E = add(scale(mul(mul(power(km, 3), z), M), 27),
            scale(mul(mul(power(k, 2), power(km, 2)), D), 540), scale(power(D, 4), -15))
    Tbar = add(mul(k, add(scale(mul(power(k, 2), power(z, 2)), -180),
                         scale(mul(D, power(z, 4)), 30), scale(mul(km, power(z, 5)), 18))),
               scale(power(k, 5), 60), scale(mul(D, power(k, 5)), -6),
               scale(mul(km, power(k, 6)), 3))
    G = add(mul(add(scale(mul(k, z), 6), scale(power(k, 2), 3)), mul(z, M)), Tbar)
    B = {}
    for (i, j), c in Q.items():
        B[(i, 0)] = B.get((i, 0), F(0)) + 64*c*F(9, 4)**j
    result = {'N': N, 'E': E, 'G': G, 'upper_z': B}
    for name, p in result.items():
        assert all(F(c).denominator == 1 for c in p.values()), name
        result[name] = {m: int(c) for m, c in p.items() if c}
    return result


def identity_check():
    # Verify the elimination of P'(-1) at the first minimum in an
    # independent polynomial ring with variables g,z,u,k.
    from hermite_criterion import add as a, mul as m, scale as sc, power as pw
    one = {(0, 0, 0, 0): F(1)}
    g, z, u, k = ({(1, 0, 0, 0): F(1)}, {(0, 1, 0, 0): F(1)},
                  {(0, 0, 1, 0): F(1)}, {(0, 0, 0, 1): F(1)})
    power3 = lambda p, n: pw(p, n, 4)
    h = m(g, z)
    alpha = m(power3(g, 3), u)
    beta = m(power3(g, 2), a(sc(one, 3), g))
    qminus = a(m(alpha, h), sc(m(beta, power3(h, 2)), F(-1, 2)),
               sc(power3(h, 4), F(1, 4)), sc(power3(h, 5), F(1, 20)))
    difference = a(sc(m(qminus, h), -1), sc(m(alpha, power3(h, 2)), F(1, 2)),
                   sc(m(beta, power3(h, 3)), F(-1, 6)), sc(power3(h, 5), F(1, 20)),
                   sc(power3(h, 6), F(1, 120)))
    J = a(sc(m(z, a(sc(one, 3), g)), F(1, 3)), sc(power3(z, 3), F(-1, 5)),
          sc(m(g, power3(z, 4)), F(-1, 24)))
    assert difference == m(m(power3(g, 5), power3(z, 2)), a(J, sc(u, F(-1, 2))))
    f = m(k, g)
    difference_c = a(m(qminus, f), sc(m(alpha, power3(f, 2)), F(1, 2)),
                     sc(m(beta, power3(f, 3)), F(1, 6)), sc(power3(f, 5), F(-1, 20)),
                     sc(power3(f, 6), F(1, 120)))
    H = a(m(k, z), sc(power3(k, 2), F(1, 2)))
    T = a(m(k, a(sc(m(a(sc(one, 3), g), power3(z, 2)), F(-1, 2)),
                     sc(power3(z, 4), F(1, 4)), sc(m(g, power3(z, 5)), F(1, 20)))),
          sc(m(a(sc(one, 3), g), power3(k, 3)), F(1, 6)),
          sc(power3(k, 5), F(-1, 20)), sc(m(g, power3(k, 6)), F(1, 120)))
    assert difference_c == m(power3(g, 5), a(m(H, u), T))


def encode(p):
    return [[list(m), c] for m, c in sorted(p.items())]


def make_certificate():
    ps = polynomials()
    return {'mask': MASK, 'polynomials': {n: encode(p) for n, p in ps.items()},
            'boxes': [{'k': [[1, 1], [11, 8]], 'z': [[1, 1], [2, 1]], 'polynomial': 'G'},
                      {'k': [[1, 1], [19, 16]], 'z': [[2, 1], [9, 4]], 'polynomial': 'N'},
                      {'k': [[19, 16], [41, 32]], 'z': [[2, 1], [9, 4]], 'polynomial': 'N'},
                      {'k': [[41, 32], [11, 8]], 'z': [[2, 1], [9, 4]], 'polynomial': 'E'}]}


def check_certificates():
    identity_check()
    expected = make_certificate()
    data = json.loads((ROOT/'data/quintic_minimum_certificate.json').read_text())
    assert data == expected
    ps = polynomials()
    bb = bernstein2(dense(ps['upper_z']), F(1), F(11, 8), F(2), F(9, 4))
    assert min(c for row in bb for c in row) > 0
    minima = []
    for box in data['boxes']:
        lo, hi = [F(*v) for v in box['k']]
        zl, zh = [F(*v) for v in box['z']]
        bs = bernstein2(dense(ps[box['polynomial']]), lo, hi, zl, zh)
        low = min(c for row in bs for c in row)
        assert low >= 0
        minima.append(str(low))
    print('EXACTLY VERIFIED two value identities, upper z bound, and 4 Bernstein rectangles;',
          'coefficient minima:', minima, flush=True)
    return data


def reasons(word):
    sig = patterns(word)[0][0]
    return ['quintic_first_minimum'] if all(sig[i] == v for i, v in MASK) else []


if __name__ == '__main__':
    import sys
    if '--write' in sys.argv:
        (ROOT/'data/quintic_minimum_certificate.json').write_text(json.dumps(make_certificate(), indent=2)+'\n')
    check_certificates()
