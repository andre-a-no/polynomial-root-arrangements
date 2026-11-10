"""Command-line entry point; no third-party Python dependency is needed."""

import argparse
import json
from pathlib import Path
from .backend import ROOT
from .pipeline import verify_all
from .stages import STAGES
from .witnesses import read_json, verify_collection, verify_witness
from .words import count_words, enumerate_words


def main():
    parser = argparse.ArgumentParser(description="Enumerate and exactly verify degree-six Rolle words")
    commands = parser.add_subparsers(dest="command", required=True)
    enumerate_parser = commands.add_parser("enumerate", help="Enumerate formal Rolle words")
    enumerate_parser.add_argument("--degree", type=int, default=6)
    enumerate_parser.add_argument("--output", type=Path)
    verify_parser = commands.add_parser("verify", help="Run every proof and realization check")
    verify_parser.add_argument("--output-dir", type=Path, default=Path("verification_output"))
    commands.add_parser("stages", help="List the proof stages, without claiming verification")
    witness_parser = commands.add_parser("witness", help="Verify and display one stored realization")
    witness_parser.add_argument("word")
    witness_parser.add_argument("--degree", type=int, choices=(5, 6), default=6)
    commands.add_parser("verify-witnesses", help="Check realizations only; this is a lower-bound check")
    args = parser.parse_args()
    if args.command == "enumerate":
        print(f"Degree {args.degree}: {count_words(args.degree)} formal words")
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(list(enumerate_words(args.degree)), indent=2) + "\n")
    elif args.command == "verify":
        verify_all(args.output_dir)
    elif args.command == "stages":
        for number, stage in enumerate(STAGES, start=1):
            print(f"{number:02} {stage.expected_new:4}  {stage.title} [{stage.module}]")
        print("Counts above are reference expectations, not a verification result.")
    elif args.command == "witness":
        witnesses = read_json(ROOT / f"data/certificates_n{args.degree}.json")
        if args.word not in witnesses:
            parser.error("No realization certificate is stored for that word")
        verify_witness(args.word, witnesses[args.word], args.degree)
        print(json.dumps({"status": "verified", "word": args.word,
                          "certificate": witnesses[args.word]}, indent=2))
    else:
        for degree in (5, 6):
            words = verify_collection(ROOT / f"data/certificates_n{degree}.json", degree)
            print(f"Degree {degree}: {len(words)} exact realizations verified")


if __name__ == "__main__":
    main()
