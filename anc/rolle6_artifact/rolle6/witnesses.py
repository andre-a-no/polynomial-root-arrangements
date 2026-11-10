"""Small integer-only verifier for realizability certificates.

A degree-d polynomial with opposite signs on d disjoint open intervals
has one simple root in each interval and no other root. Applying this
to all derivatives proves both hyperbolicity and the complete root order.
"""

import json
from pathlib import Path
from .words import validate_word


class VerificationError(ValueError):
    """A mathematical or bookkeeping verification failed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_unique_object)


def from_integer_roots(roots: list[int]) -> list[int]:
    """Expand product(x-r), using descending coefficient order."""
    coefficients = [1]
    for root in roots:
        expanded = [0] * (len(coefficients) + 1)
        for index, coefficient in enumerate(coefficients):
            expanded[index] += coefficient
            expanded[index + 1] -= root * coefficient
        coefficients = expanded
    return coefficients


def differentiate(coefficients: list[int]) -> list[int]:
    degree = len(coefficients) - 1
    return [(degree - index) * value for index, value in enumerate(coefficients[:-1])]


def scaled_value(coefficients: list[int], numerator: int, denominator: int) -> int:
    """Evaluate denominator**degree * p(numerator/denominator) exactly."""
    value = coefficients[0]
    power = denominator
    for coefficient in coefficients[1:]:
        value = value * numerator + coefficient * power
        power *= denominator
    return value


def verify_witness(word: str, item: dict, degree: int) -> None:
    validate_word(word, degree)
    require(isinstance(item, dict), "Witness must be an object")
    if "roots" in item:
        roots = item["roots"]
        require(isinstance(roots, list) and len(roots) == degree, "Invalid root list")
        require(all(type(root) is int for root in roots), "Roots must be integers")
        require(all(a < b for a, b in zip(roots, roots[1:])), "Roots must be distinct")
        coefficients = from_integer_roots(roots)
    else:
        coefficients = item.get("coefficients")
        require(isinstance(coefficients, list) and len(coefficients) == degree + 1,
                "Invalid coefficient list")
        require(all(type(value) is int for value in coefficients), "Noninteger coefficient")
        require(coefficients[0] != 0, "Incorrect degree")

    denominator = item.get("denominator")
    boundaries = item.get("boundaries")
    require(type(denominator) is int and denominator > 0, "Invalid denominator")
    require(isinstance(boundaries, list) and len(boundaries) == len(word) + 1,
            "Invalid interval boundaries")
    require(all(type(value) is int for value in boundaries), "Noninteger boundary")
    require(all(a < b for a, b in zip(boundaries, boundaries[1:])),
            "Intervals must have disjoint interiors and ordered boundaries")
    derivatives = [coefficients]
    for _ in range(degree - 1):
        derivatives.append(differentiate(derivatives[-1]))
    for index, label in enumerate(word):
        polynomial = derivatives[int(label)]
        left = scaled_value(polynomial, boundaries[index], denominator)
        right = scaled_value(polynomial, boundaries[index + 1], denominator)
        require(left * right < 0, f"{word}: interval {index} has no certified sign change")


def verify_collection(path: Path, degree: int) -> set[str]:
    data = read_json(path)
    require(isinstance(data, dict), "Witness collection must be an object")
    for word, item in data.items():
        verify_witness(word, item, degree)
    return set(data)
