"""Direct sign consequences of exact endpoint quadrature identities.

See PROOF_MARKER_QUADRATURE.md. Reflection and strict marker signs use
the same canonical projection as the previously verified obstructions.
This is a separate stage to preserve the historical exclusion counts.
"""
from inflection_obstruction import patterns
from laguerre_trapezoid_obstruction import check_identities

CONDITIONS = {
    'outer_cubic_crossing': [(0, 1), (1, 1), (15, -1), (16, 1)],
    'fourth_derivative_crossing': [(3, 1), (4, 1), (12, -1), (13, 1)],
}

def reasons(word):
    sig = patterns(word)[0][0]
    return ['marker_quadrature_' + name for name, mask in CONDITIONS.items()
            if all(sig[i] == sign for i, sign in mask)]

if __name__ == '__main__':
    check_identities()
