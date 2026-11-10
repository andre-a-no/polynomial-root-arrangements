# A readable enumeration and verification pipeline

## Mathematical objects

A degree-n word uses digit j exactly n-j times. For every prefix, if c[j]
counts digit j, adjacent interlacing is equivalent to

```text
0 <= c[j] - c[j+1] <= 1  for every j < n-1.
```

`words.enumerate_words` appends a digit only when these inequalities and the
target multiplicities remain valid. Each word is yielded once, in lexicographic
order. `words.count_words` memoizes counts of completions instead of storing
words. Labels are single digits; the interface accepts degrees 1 through 9.
Enumerating all words in the higher degrees can be expensive.

## Realizations

`data/certificates_n6.json` is an object keyed by the 2706 distinct words.
Each value gives an integer polynomial, a positive integer `denominator`, and
22 increasing integer `boundaries`. Polynomial coefficients are in descending
order; an integer `roots` list is an alternative representation of a monic
polynomial. If both representations are present, `roots` is authoritative.

For letter j at position i, the verifier checks opposite signs of P^(j) at
boundaries[i]/denominator and boundaries[i+1]/denominator. The degree and the
number of disjoint sign-changing intervals prove all roots and their order.
`scaled_value` evaluates q^d f(b/q) by integer Horner arithmetic. There are no
root-finding tolerances. Duplicate JSON object keys are rejected.

## Exclusions

Every entry in `stages.STAGES` names:

- a stable key and readable title;
- an exact certificate checker in the archived backend;
- a predicate, or a way to interpret the checker's returned words;
- the mathematical source note;
- the expected incremental count and, when available, the archived word list.

The expected count never selects a word. For each rule the pipeline first
executes its proof checker. It then computes all matches among the original
13404 candidates and checks that none is a verified realization. Only matches
still present are removed. Each new exclusion is assigned exactly one rule.

`verification.json` contains `new_exclusions` for each stage, `excluded_by_rule`
for every excluded word, and `realizable_words` for every survivor. It records
overlap as `already_excluded`. The final check compares full sets, not counts.

The initial necessary filters are distinct from the 53 stages: delete all
zeros, lower the other labels by one, and require one of the 116 certified
quintic words; then apply Kostov--Shapiro to each applicable derivative and
reflection. The published theorem that the quintic count is 116 supplies
completeness of that stored list.

## Trust boundary and preservation

`backend.py` imports the specialized modules only from this package's
`research/` directory. Their exact certificate checks are run before exclusions;
identical checker calls are cached only within that run. Search routines are
not called. Some retained modules also contain historical search entry points;
those are outside the verification path.

The mathematical modules and data are unchanged from the completed research
archive. `release/legacy_hashes.json` records their frozen SHA-256 hashes.
`dependencies.json` records the inputs actually used. Hashes detect changes;
they do not prove mathematical soundness. The archived backend uses assertions,
so `verify` refuses optimized Python (`-O` / `PYTHONOPTIMIZE`). New front-end
checks use explicit exceptions.

The compact public artifact contains the dependencies of a successful full
run, the proof notes referenced by the stage registry, and the earlier
enumerator needed for equivalence tests. Historical failed searches and
polynomial caches are not needed.

## Meaningful checks

The tests compare every formal word and every published-filter decision with
the previous implementation for degrees 1-6; check integer evaluation against
rational arithmetic; verify all 2822 quintic and sextic witnesses; reject
corrupted certificates and duplicate keys; and ensure optimized execution
cannot bypass assertions. The full integration run verifies all 53 stages.

These tests address implementation and packaging errors. They cannot replace
human assessment of the necessary inequalities, published assumptions, or
novelty of the result.
