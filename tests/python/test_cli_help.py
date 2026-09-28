"""``--help`` is documentation, and it is the copy a user reads before any other.

Two things it has to get right, both of which it got wrong before these tests existed: the column
names in it must render verbatim (rich rewrites ``:mask:`` to a pictograph, so
``vsig:mask:TRB:present`` came out as ``vsig<emoji>*``), and the corpus list must be the one the
library actually ships (it named 2 of 9 for a whole release).
"""
from __future__ import annotations

from typer.testing import CliRunner

from vdjtools.cli import app
from vdjtools.signature.corpus import bundled_names

runner = CliRunner()


def _help(*args: str) -> str:
    result = runner.invoke(app, [*args, "--help"], env={"COLUMNS": "200"})
    assert result.exit_code == 0, result.output
    return result.output


def test_the_signature_help_prints_column_names_verbatim_not_as_emoji():
    out = _help("signature")
    assert "vsig:mask:" in out, "the mask channel prefix must survive rich's emoji substitution"
    # The box-drawing frame is legitimately non-ASCII; an emoji is not. Pin the one that bit us.
    assert "\U0001f637" not in out


def test_the_corpus_flag_help_names_every_published_corpus():
    out = " ".join(_help("signature").split())
    missing = [n for n in bundled_names() if n not in out]
    assert not missing, f"--corpus help does not name published corpora: {sorted(missing)}"


def _documented_commands() -> set[str]:
    """Command names appearing in ``docs/cli.rst``, from its literal blocks and tables."""
    import re
    from pathlib import Path

    doc = Path(__file__).resolve().parents[2] / "docs" / "cli.rst"
    text = doc.read_text()
    # ``vdjtools <cmd>`` in a shell example, or ``cmd`` / ``cmd`` / ``cmd`` in a heading or table.
    return set(re.findall(r"vdjtools (\w[\w-]*)", text)) | set(re.findall(r"``([\w-]+)``", text))


def test_every_command_is_in_the_command_reference():
    """A reference page that silently omits a command is worse than no page.

    ``tcrnet``, ``alice``, ``model net``, ``model entropy``, ``model build`` and ``model rescale``
    were all shipped and undocumented before this test existed.
    """
    from typer.main import get_command

    documented = _documented_commands()
    root = get_command(app)
    missing = []
    for name, sub in root.commands.items():
        inner = getattr(sub, "commands", None)
        if inner:
            missing += [f"{name} {s}" for s in inner if s not in documented]
        elif name not in documented:
            missing.append(name)
    assert not missing, f"docs/cli.rst does not mention: {sorted(missing)}"
