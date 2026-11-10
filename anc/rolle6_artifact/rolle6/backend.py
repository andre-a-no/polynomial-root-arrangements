"""Bridge to the archived, unchanged exact mathematical verifiers.

Search programs are never imported. Proof checks run before a rule may
remove any word. Identical checker calls are memoized within one run only.
"""

from contextlib import redirect_stdout
from importlib import import_module
from pathlib import Path
import sys
from .witnesses import require

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research"


class ProofBackend:
    def __init__(self, log):
        require(__debug__, "Certificate verification requires assertions: do not use python -O")
        self.log = log
        self.cache = {}
        sys.path.insert(0, str(RESEARCH))

    def module(self, name):
        module = import_module(name)
        require(Path(module.__file__).resolve().parent == RESEARCH,
                f"Unexpected proof-module location: {name}")
        return module

    def checked(self, module_name, function_name, arguments=()):
        key = (module_name, function_name, arguments)
        if key not in self.cache:
            print(f"\n--- {module_name}.{function_name}{arguments} ---", file=self.log, flush=True)
            with redirect_stdout(self.log):
                self.cache[key] = getattr(self.module(module_name), function_name)(*arguments)
            self.log.flush()
        return self.cache[key]
