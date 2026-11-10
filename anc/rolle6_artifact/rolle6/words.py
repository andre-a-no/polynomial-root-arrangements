"""Rolle words as interlacing prefixes, without numerical root finding.

If c[j] counts occurrences of j in a prefix, adjacent derivative chains
interlace exactly when 0 <= c[j] - c[j+1] <= 1 for every prefix.
The final multiplicities are (n, n-1, ..., 1).
"""

from functools import lru_cache
from typing import Iterator


def _multiplicities(degree: int) -> tuple[int, ...]:
    if type(degree) is not int or not 1 <= degree <= 9:
        raise ValueError("Single-digit Rolle words require 1 <= degree <= 9")
    return tuple(range(degree, 0, -1))


def _successors(counts: tuple[int, ...], limits: tuple[int, ...]):
    for label, count in enumerate(counts):
        if count == limits[label]:
            continue
        # Adding label must not overtake its predecessor or leave its
        # successor more than one root behind.
        if label > 0 and count == counts[label - 1]:
            continue
        if label + 1 < len(counts) and count > counts[label + 1]:
            continue
        updated = list(counts)
        updated[label] += 1
        yield label, tuple(updated)


def enumerate_words(degree: int) -> Iterator[str]:
    """Yield every admissible word once, in lexicographic order."""
    limits = _multiplicities(degree)

    def visit(counts: tuple[int, ...], prefix: str) -> Iterator[str]:
        if counts == limits:
            yield prefix
            return
        for label, updated in _successors(counts, limits):
            yield from visit(updated, prefix + str(label))

    yield from visit((0,) * degree, "")


def count_words(degree: int) -> int:
    """Count extensions by dynamic programming, without storing the words."""
    limits = _multiplicities(degree)

    @lru_cache(maxsize=None)
    def count(counts: tuple[int, ...]) -> int:
        if counts == limits:
            return 1
        return sum(count(updated) for _, updated in _successors(counts, limits))

    return count((0,) * degree)


def validate_word(word: str, degree: int) -> None:
    limits = _multiplicities(degree)
    if not isinstance(word, str) or len(word) != sum(limits):
        raise ValueError("Incorrect word length")
    counts = [0] * degree
    for character in word:
        if character not in "012345678"[:degree]:
            raise ValueError("Invalid derivative label")
        label = int(character)
        counts[label] += 1
        if counts[label] > limits[label]:
            raise ValueError("Incorrect derivative multiplicity")
        if any(not 0 <= counts[j] - counts[j + 1] <= 1 for j in range(degree - 1)):
            raise ValueError("Prefix violates strict interlacing")
    if tuple(counts) != limits:
        raise ValueError("Incorrect final multiplicities")


def derivative_word(word: str) -> str:
    return "".join(str(int(label) - 1) for label in word if label != "0")


def root_positions(word: str, degree: int) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(i for i, label in enumerate(word) if label == str(j))
                 for j in range(degree))


def kostov_shapiro_admissible(word: str, degree: int) -> bool:
    """Apply Theorem 1 of Kostov--Shapiro to each derivative and reflection.

    If every internal root of f precedes the matching root of f'', the
    middle root of f^(m-3) must precede the root of f^(m-1). Reflection
    gives the reversed implication. Indices here start at zero.
    """
    positions = root_positions(word, degree)
    for derivative in range(degree - 3):
        size = degree - derivative
        comparisons = [positions[derivative][j + 1] < positions[derivative + 2][j]
                       for j in range(size - 2)]
        middle_precedes_mean = positions[degree - 3][1] < positions[degree - 1][0]
        if all(comparisons) and not middle_precedes_mean:
            return False
        if not any(comparisons) and middle_precedes_mean:
            return False
    return True


def primitive_mean_signature(word: str) -> str:
    # For a sextic with positive leading coefficient, the sign at its
    # root barycentre is (-1)^(number of roots to its left).
    centre = word.index("5")
    sign = (-1) ** word[:centre].count("0")
    return f"{derivative_word(word)}:{sign}"
