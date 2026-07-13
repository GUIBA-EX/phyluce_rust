"""Locate optional legacy sources and checked-in compatibility fixtures."""

import os
from pathlib import Path

RUST_ROOT = Path(__file__).resolve().parents[2]
LOCAL_FIXTURE_REPO = RUST_ROOT / "tests/compat/fixtures/python-repo"


def _is_python_repo(path: Path) -> bool:
    return (path / "bin").is_dir() and (path / "phyluce").is_dir()


def find_python_repo() -> Path:
    configured = os.environ.get("PHYLUCE_PYTHON_REPO")
    if configured:
        path = Path(configured).expanduser().resolve()
        if _is_python_repo(path):
            return path
        raise RuntimeError(f"PHYLUCE_PYTHON_REPO is not a phyluce source tree: {path}")

    nested_checkout = RUST_ROOT.parent
    if _is_python_repo(nested_checkout):
        return nested_checkout
    raise RuntimeError(
        "live Python/Rust comparison requires the original phyluce source; "
        "set PHYLUCE_PYTHON_REPO=/path/to/phyluce"
    )


def find_fixture_repo() -> Path:
    if LOCAL_FIXTURE_REPO.is_dir():
        return LOCAL_FIXTURE_REPO
    return find_python_repo()
