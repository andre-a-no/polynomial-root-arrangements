# Prepared arXiv source package

Title: An exact certificate-based enumeration of Rolle words of degree six

Author: Andrey A. Novikov

Version: 1.0, 10 November 2026

Suggested primary category: math.CA. Suggested cross-list: math.CO.
These are preparation suggestions, not a record of a submission.

Suggested comments: 11 pages; exact certificate corpus and Python verification
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

For a strictly hyperbolic real polynomial whose nonconstant derivatives have
mutually distinct roots, record the derivative orders of all roots from left
to right. We report a computer-assisted classification of these Rolle words
in degree six: exactly 2706 words are realizable, counting reflections
separately. The computation starts from 33592 formal interlacing words. The
known quintic classification and a theorem of Kostov--Shapiro reduce this to
13404 candidates. Integer realization certificates and 53 exact exclusion
stages produce an exhaustive, disjoint classification. Verification uses
integer and rational arithmetic; floating-point searches are not used as
impossibility proofs. A readable enumerator and an explicit proof-stage ledger
accompany the certificate corpus. Most of the mathematical development,
software, certificate construction, and initial exposition were generated
using ChatGPT in an interactive project directed by the author. This
substantial contribution and the absence of an independent human audit are
documented explicitly.
