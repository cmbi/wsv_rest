from wiws.interface import run_whatif


def test_whatif():
    lines = run_whatif("WSVDOT", "2WDQ", False)
    assert len(lines) > 0
