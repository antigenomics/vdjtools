"""Every release tag must have a CHANGELOG section, and every section a tag.

This is a test because the drift is silent in both directions and both directions happened.
**v3.9.2 and v3.17.0 were tagged, built and published with no section at all** -- found in 4.5.0 by
reading the tag list against the file, not by anyone noticing. A release whose only record is a
commit message is a release nobody can read the history of, and the gap is invisible from inside the
file: the headings above and below it look perfectly ordinary.

The other direction catches a typo in a heading, which presents as a missing section for the real
version plus a section for a version that does not exist.

Deliberately not asserted: that a section's date equals its tag's commit date. Releases land across
midnight -- `v3.15.0`'s commits are 00:01 and 00:13 -- so the section's date is the day the work was
done and a strict comparison would fail on a true record.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]

#: v2 releases before this are recorded in the git tags only, as the file's own preamble says.
FIRST_DOCUMENTED = (3, 0, 0)

SECTION = re.compile(r"^## (\d+\.\d+\.\d+)(?:\s+—\s+(\d{4}-\d{2}-\d{2}))?\s*$", re.M)


def _parts(v: str) -> tuple[int, ...]:
    return tuple(int(x) for x in v.split("."))


@pytest.fixture(scope="module")
def versions() -> tuple[dict[str, str], set[str]]:
    """Sections mapped to their date, and the tagged versions this file is expected to cover."""
    r = subprocess.run(["git", "tag", "--list", "v*.*.*"], cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        pytest.skip("not a git checkout, so there are no tags to compare against")
    sections = dict(SECTION.findall((ROOT / "CHANGELOG.md").read_text()))
    tags = {t[1:] for t in r.stdout.split() if re.fullmatch(r"v\d+\.\d+\.\d+", t)}
    assert tags, "no release tags found; the comparison would pass by being empty"
    return sections, {t for t in tags if _parts(t) >= FIRST_DOCUMENTED}


def test_every_release_tag_has_a_section(versions):
    sections, tags = versions
    missing = sorted(tags - set(sections), key=_parts)
    assert not missing, (
        f"tagged and released with no CHANGELOG section: {missing}. Write it from the release "
        f"commit (`git log v<prev>..v<ver>`), not from memory."
    )


def test_every_section_names_a_real_tag(versions):
    """Every section but the newest, which may be a release in flight.

    The exemption is not a loophole, it is the workflow: a release writes its section, bumps the
    version and commits, and only then tags. Without it this test goes red on every release
    between those two steps -- which it did, on 4.6.0, the first release after it was written.
    Only the single highest version is exempt, so a typo'd heading is still caught the moment
    anything lands above it, and a *stale* untagged section is caught by the next release.
    """
    sections, tags = versions
    unknown = set(sections) - tags
    if unknown:
        newest = max(sections, key=_parts)
        unknown -= {newest}                      # the release being prepared right now
    assert not unknown, (
        f"CHANGELOG sections for versions that were never tagged: {sorted(unknown, key=_parts)}. "
        f"Usually a typo in the heading, which also hides the real version as a missing section."
    )


def test_every_release_section_is_dated(versions):
    sections, _ = versions
    undated = sorted((v for v, date in sections.items() if not date), key=_parts)
    assert not undated, f"release sections with no date: {undated}"


def test_there_is_at_most_one_unreleased_section_and_it_comes_first(versions):
    """A second one buried mid-file is work that was released without being attributed to anything.

    That is what happened to the seqtree floor note: it sat between the 3.14.0 and 3.13.0 entries
    for five releases while the change itself had shipped in 3.14.0.
    """
    text = (ROOT / "CHANGELOG.md").read_text()
    heads = re.findall(r"^## (.+)$", text, re.M)
    unreleased = [i for i, h in enumerate(heads) if h.strip().lower() == "unreleased"]
    assert len(unreleased) <= 1, (
        f"{len(unreleased)} '## Unreleased' sections. A later one is released work with no version "
        f"against it -- fold it into the release that shipped it."
    )
    if unreleased:
        assert unreleased[0] == 0, f"'## Unreleased' must be the first section, found after {heads[0]!r}"
