# Verification record — 10 November 2026

The complete refactored verification passed in the research tree and again
inside this compact ancillary artifact. The two full per-word stage ledgers
are identical, after ignoring elapsed times.

| Check | Result |
| --- | ---: |
| Formal degree-six words | 33592 |
| After the complete quintic-tail filter | 13950 |
| After the Kostov–Shapiro filter | 13404 |
| Exact exclusion stages | 53 |
| Distinct excluded candidates | 10698 |
| Exactly verified sextic realizations | 2706 |
| Unresolved candidates | 0 |
| Exactly verified quintic realizations | 116 |
| Refactoring tests | 6 passed |
| First full refactored run | 483.665 s |
| Full compact-artifact run | 480.025 s |

Python used in this preparation environment: 3.12.14.
Timing is descriptive and is not a performance guarantee.

Files:

- `release/clean_verification/verification.json`: complete per-word accounting.
- `release/clean_verification/dependencies.json`: actual input hashes.
- `release/clean_verification/proof_checkers.log`: specialized checker outputs.
- `release/clean_verification.log`: high-level execution log.
- `release/unit_tests.log`: equivalence and damaged-certificate test results.
- `release/legacy_hashes.json`: frozen mathematical sources and data.

The final surviving **set** equals the verified witness set. Every new
exclusion is disjoint from earlier new exclusions. The mathematical proof
modules and data are unchanged from the completed research archive.

No independent human audit of the entire mathematical proof chain has been
completed. The new front end and the historical implementation share the
specialized proof backend. Published dependencies remain Kostov's quintic
classification and the Kostov–Shapiro restriction. Tests and exact execution
do not by themselves establish novelty.
