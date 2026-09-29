"""Working notes must never enter the public repository or a published artifact.

Several markdown files in this tree exist to hold internal detail -- cluster project paths,
private HuggingFace dataset names, cohort composition, in-flight decisions. ``SOURCES.md`` is the
clearest case: its *purpose* is provenance that must not be published. ``CLAUDE.md``, ``ISSUES.md``
and ``TODO.md`` are the same in kind.

They are kept out by two independent mechanisms, and this pins both:

* ``.gitignore`` keeps them untracked, so ``git add -A`` cannot sweep one in;
* ``pyproject.toml``'s sdist exclude keeps them out of the published package.

Neither is self-announcing when it breaks -- a newly tracked ``NOTES.md`` looks like any other
file in a diff -- which is why this is a test rather than a convention.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

try:
    import tomllib                    # stdlib from 3.11
except ModuleNotFoundError:           # 3.10: declared in the [test] extra
    import tomli as tomllib

import pytest

ROOT = Path(__file__).resolve().parents[2]

# Filenames whose content is internal by design. Matched against the basename, case-insensitively.
PRIVATE = re.compile(
    r"^(SOURCES|CLAUDE|TODO|NOTES|ROADMAP|STATUS|AGENTS|NULLS|PLAN)\.md$|^ISSUES.*\.md$",
    re.IGNORECASE,
)


# Every working-note name, as the sdist exclude and the gitignore must both spell them.
NAMES = ("SOURCES.md", "CLAUDE.md", "TODO.md", "NOTES.md", "ROADMAP.md", "STATUS.md",
         "ISSUES.md", "ISSUES_ext.md")


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout


def test_no_working_note_is_tracked_by_git():
    tracked = [p for p in _git("ls-files").split("\n") if p]
    leaked = sorted(p for p in tracked if PRIVATE.match(Path(p).name))
    assert not leaked, (
        f"working notes are tracked and would be published: {leaked}. "
        f"Run `git rm --cached <file>` and confirm .gitignore covers it."
    )


@pytest.mark.parametrize("name", NAMES)
def test_each_working_note_name_is_gitignored(name):
    """Ignored by NAME, not merely absent -- an absent file that is not ignored is a future leak."""
    r = subprocess.run(["git", "check-ignore", "-q", name], cwd=ROOT)
    assert r.returncode == 0, (
        f"{name} is not gitignored. It may not exist today, but nothing stops it being created "
        f"and committed tomorrow."
    )


def test_the_sdist_exclude_names_every_working_note():
    """Second line of defence, deliberately backend-agnostic.

    The build backend selects files from the VCS listing, so ``.gitignore`` is what actually keeps
    these out -- which is the test above. This one pins the explicit exclude as well, because that
    is the only thing standing between an untracked copy sitting in a working tree and an sdist.
    Asserted against the file's text rather than a key path: scikit-build-core and hatchling spell
    the table differently, and a test that follows one of them silently stops checking if the
    backend changes.
    """
    text = (ROOT / "pyproject.toml").read_text()
    tomllib.loads(text)          # it must still be valid TOML after any edit to that list
    missing = [n for n in NAMES if f'"{n}"' not in text and f'"/{n}"' not in text]
    assert not missing, (
        f"not excluded from the sdist: {missing}. gitignore does not reach the build backend for "
        f"a file that is present but untracked."
    )


def test_the_matcher_can_actually_fail():
    """A guard: the pattern must match the names it exists to catch, and leave public files alone."""
    assert PRIVATE.match("SOURCES.md") and PRIVATE.match("ISSUES_ext.md")
    assert PRIVATE.match("claude.md"), "matching must be case-insensitive"
    assert not PRIVATE.match("README.md")
    assert not PRIVATE.match("CHANGELOG.md")
    assert not PRIVATE.match("SKILL.md")
