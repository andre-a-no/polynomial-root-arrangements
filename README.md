# Degree-six polynomial root arrangements

**Preprint, version 1.0 (10 November 2026). Not peer reviewed.**

An exact certificate-based enumeration of Rolle words of degree six  
Human author: Andrey A. Novikov.

**Most of the mathematical development, code, exact certificates, and manuscript were generated using ChatGPT (OpenAI), in an interactive investigation directed by the human author.** ChatGPT is disclosed as a tool, not a coauthor. The human author remains responsible for any submitted version. See [the full provenance disclosure](anc/rolle6_artifact/AI_PROVENANCE.md).

- [Read the preprint (PDF)](main.pdf)
- [LaTeX source](main.tex)
- [Code, certificates, and reproduction instructions](anc/rolle6_artifact/README.md)
- [Submission preparation notes](README_submission.md)

The reported result is 2706 realizable degree-six Rolle words, with reflections counted separately. This is not a count of connected components. The exact checker supports reproducibility; no independent human audit of the complete mathematical proof chain has been completed. Publication priority has not been established.

## Reproduce

Python 3.10 or later; no third-party Python dependencies; do not use `-O`.

```sh
cd anc/rolle6_artifact
python3 -m rolle6 verify --output-dir verification_output
python3 -m unittest discover -s tests -v
```

Build the paper with pdfLaTeX, BibTeX, then pdfLaTeX twice. A standard TeX Live installation is sufficient; no shell escape is needed.

## Scope and status

This repository contains the compact preprint-and-code upload, including its required research proof modules. The separate historical research ZIP could not be inspected or imported because it exceeds the attachment transfer limit of 32 MiB. It is not represented as published here.

Prepared for a possible arXiv submission; no arXiv submission, acceptance, identifier, or endorsement is asserted. No distribution license has been selected in this release; public availability alone does not grant an open-source license.
