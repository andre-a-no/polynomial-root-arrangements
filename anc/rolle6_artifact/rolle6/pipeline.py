"""One auditable sequence: enumerate, verify witnesses, check proofs, subtract sets."""

import hashlib
import json
from pathlib import Path
import sys
import time
from .backend import ProofBackend, ROOT
from .stages import STAGES
from .witnesses import read_json, require, verify_collection
from .words import (count_words, derivative_word, enumerate_words,
                    kostov_shapiro_admissible, primitive_mean_signature)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


class ReadDependencies:
    """Record actual JSON inputs; imported source files are collected separately."""
    def __init__(self):
        self.paths = set()
        self.active = True
        sys.addaudithook(self.record)

    def record(self, event, arguments):
        if self.active and event == "open" and isinstance(arguments[0], str):
            path = Path(arguments[0])
            if path.suffix == ".json" and path.is_absolute() and ROOT in path.parents:
                if path.parent == ROOT / "data":
                    self.paths.add(path)

    def finish(self):
        self.active = False
        for module in tuple(sys.modules.values()):
            filename = getattr(module, "__file__", None)
            if filename:
                path = Path(filename).resolve()
                if ROOT in path.parents and path.suffix == ".py":
                    self.paths.add(path)
        return {str(path.relative_to(ROOT)): digest(path) for path in sorted(self.paths)}


def verify_all(output_directory: Path) -> dict:
    require(__debug__, "Do not use python -O: archived proof checkers use assertions")
    output_directory.mkdir(parents=True, exist_ok=True)
    tracker = ReadDependencies()
    started = time.monotonic()
    print("Verifying integer realizability certificates...", flush=True)
    five = verify_collection(ROOT / "data/certificates_n5.json", 5)
    six = verify_collection(ROOT / "data/certificates_n6.json", 6)
    require(len(five) == 116 and len(six) == 2706, "Unexpected witness inventory")

    formal = set(enumerate_words(6))
    require(len(formal) == count_words(6) == 33592, "Enumeration count mismatch")
    after_tail = {word for word in formal if derivative_word(word) in five}
    candidates = {word for word in after_tail if kostov_shapiro_admissible(word, 6)}
    require(len(after_tail) == 13950 and len(candidates) == 13404,
            "Published-filter counts differ from the reference computation")
    require(candidates == set(read_json(ROOT / "data/candidates_n6.json")),
            "Candidate words differ from the archived enumeration")
    require(six <= candidates, "A verified witness was rejected by a necessary filter")
    require({word[::-1] for word in six} == six, "Witness reflection audit failed")

    remaining = set(candidates)
    ledger = []
    excluded_by = {}
    with (output_directory / "proof_checkers.log").open("w", encoding="utf-8") as log:
        backend = ProofBackend(log)
        for number, stage in enumerate(STAGES, start=1):
            print(f"[{number:02}/{len(STAGES)}] Checking {stage.title}...", flush=True)
            begin = time.monotonic()
            certificate_result = backend.checked(stage.module, stage.checker, stage.arguments)
            if stage.selection == "returned_words":
                matched = set(certificate_result)
            elif stage.selection == "tail":
                matched = {word for word in candidates if derivative_word(word) in certificate_result}
            elif stage.selection == "mean":
                matched = {word for word in candidates
                           if primitive_mean_signature(word) in certificate_result}
            else:
                predicate = getattr(backend.module(stage.module), stage.predicate)
                matched = {word for word in candidates if predicate(word)}
            require(matched <= candidates, f"{stage.key}: rule returned a noncandidate word")
            require(not matched & six, f"{stage.key}: rule contradicts a verified realization")
            new = remaining & matched
            require(len(new) == stage.expected_new, f"{stage.key}: incremental count changed")
            if stage.reference_file is not None:
                require(new == set(read_json(ROOT / stage.reference_file)),
                        f"{stage.key}: exact exclusion list changed")
            require(new == {word[::-1] for word in new}, f"{stage.key}: reflection audit failed")
            for word in new:
                require(word not in excluded_by, "Double-counted exclusion")
                excluded_by[word] = stage.key
            before = len(remaining)
            remaining -= new
            row = {"number": number, "key": stage.key, "title": stage.title,
                   "checker": f"{stage.module}.{stage.checker}", "proof_note": stage.proof_note,
                   "matched_words": len(matched), "already_excluded": len(matched - new),
                   "new_exclusions": sorted(new), "new_count": len(new),
                   "before": before, "remaining": len(remaining),
                   "seconds": round(time.monotonic() - begin, 3)}
            ledger.append(row)
            write_json(output_directory / "progress.json", {"status": "in_progress", "stages": ledger})
            print(f"           new {len(new):4}; remaining {len(remaining):5}", flush=True)

    require(remaining == six, "The final set is not exactly the set of verified witnesses")
    require(set(excluded_by) | six == candidates and not set(excluded_by) & six,
            "Final partition is not disjoint and exhaustive")
    dependencies = tracker.finish()
    frozen = read_json(ROOT / "release/legacy_hashes.json")
    for filename, checksum in dependencies.items():
        if filename.startswith(("research/", "data/")):
            require(frozen.get(filename) == checksum, f"Archived proof input changed: {filename}")

    report = {"status": "verified_exact_count", "degree": 6, "exact_count": len(six),
              "formal_words": len(formal), "after_quintic_filter": len(after_tail),
              "after_kostov_shapiro_filter": len(candidates), "excluded_words": len(excluded_by),
              "unresolved_words": 0, "stages": ledger,
              "realizable_words": sorted(six), "excluded_by_rule": dict(sorted(excluded_by.items())),
              "legacy_proof_inputs_unchanged": True,
              "external_theorems": ["Kostov (2002), Theorem 2: 116 quintic words",
                                    "Kostov--Shapiro (2002), Theorem 1"],
              "arithmetic": "Python integers and Fraction; no floating-point proof decisions",
              "seconds": round(time.monotonic() - started, 3)}
    write_json(output_directory / "verification.json", report)
    write_json(output_directory / "dependencies.json", dependencies)
    write_json(output_directory / "progress.json", {"status": "completed", "exact_count": len(six)})
    print(f"VERIFIED: a6 = {len(six)}; unresolved = 0; all {len(STAGES)} stages passed.", flush=True)
    return report
