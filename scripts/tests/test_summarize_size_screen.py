import importlib.util
from fractions import Fraction
from pathlib import Path


SCRIPT = (
    Path(__file__).parents[2]
    / ".agents/skills/boundgsenergy-hpc-operations/scripts/summarize_size_screen.py"
)
SPEC = importlib.util.spec_from_file_location("summarize_size_screen", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
SUMMARY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUMMARY)


def test_exact_replay_endpoint_prefers_rational_value():
    result = {
        "hybrid_replay": {
            "strict_lower_endpoint": 1.0,
            "strict_lower_endpoint_rational": "10000000000000001//10000000000000000",
        }
    }

    assert SUMMARY.exact_replay_endpoint(result, "hybrid_replay") == Fraction(
        10000000000000001,
        10000000000000000,
    )
