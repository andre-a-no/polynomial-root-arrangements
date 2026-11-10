# Degree-six Rolle words: exact count 2706

This is the compact reproduction artifact for the preprint. The complete
historical research archive is retained separately; this package includes only
the dependencies needed for verification, proof notes, and equivalence tests.
Historical intermediate bounds in the notes are not the final count.

The problem counts the total orders of the roots of a real degree-six
polynomial and all its nonconstant derivatives. The polynomial has six distinct
real roots; roots belonging to different derivatives must also be distinct.
Reflections are counted separately. This is not a count of connected components.

## Reproduce the result

Use **Python 3.10 or later**, without `-O`. No third-party Python dependency is
required. Run these commands from this directory:

```sh
python3 -m rolle6 enumerate --degree 6
python3 -m rolle6 verify --output-dir verification_output
python3 -m unittest discover -s tests -v
```

The first command counts 33592 formal words, not realizable words. The second
performs all witness and exclusion checks. A successful run ends with:

```text
VERIFIED: a6 = 2706; unresolved = 0; all 53 stages passed.
```

The reference run took about eight minutes in the preparation environment;
runtime depends on the machine. The output directory contains a proof-checker
log, a complete per-word exclusion ledger, the surviving witness words, and
hashes of the actual proof dependencies. A partial `progress.json` is not a
successful verification report.

Additional commands:

```sh
python3 -m rolle6 enumerate --degree 6 --output formal_words.json
python3 -m rolle6 stages
python3 -m rolle6 verify-witnesses
python3 -m rolle6 witness 010210321043210543210 --degree 6
```

For the last command, substitute any key from `data/certificates_n6.json`.
`stages` displays reference counts only. `verify-witnesses` proves the lower
bound only. Neither command substitutes for `verify`.

## Reading order

1. `rolle6/words.py`: the entire enumeration is a prefix-count recursion.
2. `rolle6/witnesses.py`: integer-only checks of explicit realizations.
3. `rolle6/stages.py`: the ordered catalogue of 53 exclusion rules.
4. `rolle6/pipeline.py`: check each proof, subtract newly excluded words, and
   compare the final set with the verified witnesses.
5. `research/`: unchanged exact mathematical checkers and their proof notes.

See `docs/ARCHITECTURE.md` for the code and certificate formats and
`AI_PROVENANCE.md` for the substantial role of ChatGPT. The preprint provides
the definitions, mathematical reductions, a worked final obstruction, and the
complete stage ledger.

## Evidence and assumptions

The computation uses Kostov's published quintic count 116 and the necessary
restriction of Kostov--Shapiro. After these filters there are 13404 candidate
words. The 53 exact stages exclude 10698 distinct words, leaving exactly the
2706 words with independently checked integer witnesses. Overlapping rules
are accounted for explicitly; no word is counted twice as a new exclusion.

The specialized proof checkers are preserved by hash from the completed
research archive. Refactoring the front end is not an independent
reimplementation of those mathematical arguments. No independent human audit
of the entire proof chain has yet been completed, and this is not a formal
proof-assistant development. Publication priority has not been established
by the computation.
