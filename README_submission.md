# Prepared arXiv source package

Title: An exact certificate-based enumeration of Rolle words of degree six

Author: Andrey A. Novikov

Version: 1.1, 11 November 2026

Suggested primary category: math.CA. Suggested cross-list: math.CO.
These are preparation suggestions, not a record of a submission.

Suggested comments: 12 pages; exact certificate corpus and Python verification
code included as ancillary files; substantial use of ChatGPT disclosed in
Section 9 and in the abstract.

## Contents

- `main.tex`: the manuscript.
- `main.bib`, `main.bbl`: bibliography source and generated bibliography.
- `stage_table.tex`: the 53-stage ledger, generated from a successful run.
- `anc/rolle6_artifact/`: self-contained code, certificates, proof notes,
  verification reports and tests.

The arXiv source ZIP contains no compiled PDF. The repository includes a separate compiled PDF for review.
The ancillary files are already in arXiv's `anc/` layout.

## Build the manuscript

In a standard TeX Live installation:

```sh
pdflatex -halt-on-error main.tex
bibtex main
pdflatex -halt-on-error main.tex
pdflatex -halt-on-error main.tex
```

The document uses standard LaTeX packages and requires no shell escape.
The supplied `main.bbl` also permits building without rerunning BibTeX.

## Reproduce the computation

```sh
cd anc/rolle6_artifact
python3 -m rolle6 verify --output-dir my_verification
python3 -m unittest discover -s tests -v
```

Use Python 3.10 or later without optimization. No third-party Python packages
are needed. See the artifact's README and architecture note for the complete
word and certificate definitions.

## Repository and release status

Public repository: https://github.com/andre-a-no/polynomial-root-arrangements

This is a preprint, not peer reviewed. The first page explicitly discloses that most of the mathematical development, code, certificates and manuscript were generated using ChatGPT under human direction.

The separate historical research ZIP was not imported: it exceeds the attachment transfer limit. The compact ancillary artifact is included in full.

## Attribution and review status

The manuscript explicitly states that ChatGPT generated most of the
mathematical development, programs, certificates and initial exposition.
Andrey A. Novikov is the human author and directed the investigation. ChatGPT
is credited as a tool, not listed as a coauthor. No independent human audit of
the entire proof chain has yet been completed. Novelty is not established by
the computation alone.

This package has not been submitted. Account details, the author's final
approval of the text, and a distribution-licence selection belong to the
actual submission step; no licence selection is asserted here.

Relevant arXiv guidance:

- https://info.arxiv.org/help/moderation/index.html
- https://info.arxiv.org/help/ancillary_files.html

## Abstract for submission

For a real polynomial of degree n with n distinct real roots, assume that no
two different nonconstant derivatives, including the polynomial itself, share
a root. Recording the derivative orders of all roots from left to right gives
a Rolle word. We present an exact, computer-assisted classification in degree
six: 2706 words are realizable, with a word and its reversal counted
separately. Starting from 33592 formal interlacing words, published necessary
restrictions leave 13404 candidates. Integer realization certificates and 53
exact exclusion stages classify every candidate. The surviving set agrees
exactly with the set of certified realizations. All verification uses integer
and rational arithmetic. The accompanying code, certificates, proof notes, and
per-word exclusion ledger make the computation reproducible. The result
concerns degree six; it does not resolve the classification problem in
arbitrary degree or count connected components of coefficient space. Most of
the mathematical development, code, certificates, and initial manuscript were
generated using ChatGPT under the author's direction. No independent human
audit of the complete proof chain has yet been completed.
