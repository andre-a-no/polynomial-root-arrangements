"""Meaningful equivalence and rejection tests for the rewritten front end."""

from copy import deepcopy
from fractions import Fraction
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from rolle6.backend import ROOT
from rolle6.witnesses import (VerificationError, read_json, scaled_value,
                              verify_collection, verify_witness)
from rolle6.words import (count_words, enumerate_words, kostov_shapiro_admissible,
                          validate_word)


class EnumerationTests(unittest.TestCase):
    def test_every_word_equals_archived_poset_enumerator(self):
        spec = importlib.util.spec_from_file_location("old_combinatorics", ROOT / "research/combinatorics.py")
        old = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(old)
        for degree, expected in enumerate((1, 1, 2, 12, 286, 33592), start=1):
            words = list(enumerate_words(degree))
            self.assertEqual(words, sorted(set(words)))
            self.assertEqual(set(words), set(old.rolle_words(degree)))
            self.assertEqual(count_words(degree), expected)
            for word in words:
                validate_word(word, degree)
                self.assertEqual(kostov_shapiro_admissible(word, degree), old.ks_ok(word, degree))

    def test_noninterlacing_word_is_rejected(self):
        with self.assertRaises(ValueError):
            validate_word("001", 2)


class WitnessTests(unittest.TestCase):
    def test_integer_evaluation_matches_rational_arithmetic(self):
        coefficients = [3, -5, 0, 11, -7]
        for denominator in (1, 3, 17):
            for numerator in (-11, -1, 0, 2, 13):
                x = Fraction(numerator, denominator)
                exact = sum(c * x ** (len(coefficients) - j - 1)
                            for j, c in enumerate(coefficients))
                self.assertEqual(scaled_value(coefficients, numerator, denominator),
                                 exact * denominator ** (len(coefficients) - 1))

    def test_all_realizations_and_corrupted_certificates(self):
        for degree, expected in ((5, 116), (6, 2706)):
            data = read_json(ROOT / f"data/certificates_n{degree}.json")
            self.assertEqual(len(verify_collection(ROOT / f"data/certificates_n{degree}.json", degree)), expected)
            word, item = next(iter(data.items()))
            for alteration in ("denominator", "boundaries", "coefficients"):
                bad = deepcopy(item)
                if alteration == "denominator":
                    bad["denominator"] = 0
                elif alteration == "boundaries":
                    bad["boundaries"][1] = bad["boundaries"][0]
                else:
                    bad.pop("roots", None)
                    bad["coefficients"] = [1] + [0] * degree
                with self.assertRaises(VerificationError):
                    verify_witness(word, bad, degree)

    def test_duplicate_json_keys_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.json"
            path.write_text('{"same": 1, "same": 2}')
            with self.assertRaises(VerificationError):
                read_json(path)

    def test_optimized_python_cannot_skip_proofs(self):
        result = subprocess.run([sys.executable, "-O", "-m", "rolle6", "verify"],
                                cwd=ROOT, text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Do not use python -O", result.stderr)


if __name__ == "__main__":
    unittest.main()
