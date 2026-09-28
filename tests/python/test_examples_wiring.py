"""Every example notebook must be internally wired before it is shipped.

marimo notebooks are plain Python: each ``@app.cell`` function declares the names it *consumes*
as parameters and the names it *provides* in its ``return`` tuple. The runtime resolves cells
against those declarations, so three mistakes are invisible to the eye and fatal at run time:

* a cell consumes a name no cell provides (an import missing from the import cell's return, or
  missing altogether);
* two cells provide the same name (marimo refuses to build the graph at all);
* a cell provides a name nothing consumes and nothing displays -- harmless, so not checked.

All three shipped in this repository at once: ``aging.py`` used ``plt`` without importing
matplotlib anywhere, ``emerson_cmv_hla.py`` used ``time`` the same way, and ``model_explorer.py``
bound ``fig`` in two cells. Each one failed on the documented ``python examples/<name>.py`` run
mode. This check is static, needs no data and no download, and runs in milliseconds.
"""

from __future__ import annotations

import ast
import builtins
import pathlib

import pytest

EXAMPLES = sorted((pathlib.Path(__file__).resolve().parents[2] / "examples").glob("*.py"))
BUILTINS = set(dir(builtins))


def _cells(tree: ast.Module) -> list[ast.FunctionDef]:
    """Every ``@app.cell``-decorated function, in file order."""
    out = []
    for node in tree.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        for dec in node.decorator_list:
            target = dec.func if isinstance(dec, ast.Call) else dec
            if isinstance(target, ast.Attribute) and target.attr == "cell":
                out.append(node)
                break
    return out


def _provided(cell: ast.FunctionDef) -> set[str]:
    """The names a cell *defines* and hands on: its ``return`` tuple minus its own parameters.

    A cell may re-return a name it received (marimo writes these pass-throughs itself), and that
    is not a second definition of it.
    """
    out: set[str] = set()
    for node in cell.body:          # top level only: a nested def has its own return
        if not isinstance(node, ast.Return) or node.value is None:
            continue
        value = node.value
        elts = value.elts if isinstance(value, (ast.Tuple, ast.List)) else [value]
        out.update(e.id for e in elts if isinstance(e, ast.Name))
    return out - _consumed(cell)


def _consumed(cell: ast.FunctionDef) -> set[str]:
    """The names a cell declares as parameters, i.e. expects another cell to provide."""
    a = cell.args
    return {p.arg for p in (*a.posonlyargs, *a.args, *a.kwonlyargs)}


assert EXAMPLES, "no example notebooks found"


@pytest.mark.parametrize("path", EXAMPLES, ids=lambda p: p.stem)
def test_every_consumed_name_is_provided_by_some_cell(path: pathlib.Path) -> None:
    cells = _cells(ast.parse(path.read_text()))
    assert cells, f"{path.name}: no @app.cell functions found"

    provided: set[str] = set()
    for cell in cells:
        provided |= _provided(cell)

    dangling = {n for cell in cells for n in _consumed(cell)} - provided - BUILTINS
    assert not dangling, (
        f"{path.name}: cell parameters that no cell returns: {sorted(dangling)}. "
        f"Either the name is never imported, or the cell that imports it left it out of its "
        f"return tuple -- both fail at run time, not at import."
    )


@pytest.mark.parametrize("path", EXAMPLES, ids=lambda p: p.stem)
def test_no_name_is_provided_by_two_cells(path: pathlib.Path) -> None:
    cells = _cells(ast.parse(path.read_text()))
    seen: dict[str, int] = {}
    clashes: dict[str, list[int]] = {}
    for i, cell in enumerate(cells):
        for name in _provided(cell):
            if name in seen:
                clashes.setdefault(name, [seen[name]]).append(i)
            else:
                seen[name] = i
    assert not clashes, (
        f"{path.name}: names returned by more than one cell: "
        f"{ {k: v for k, v in clashes.items()} }. marimo refuses to build the graph; prefix the "
        f"later one with '_' to make it cell-private."
    )


def test_the_check_can_actually_fail() -> None:
    """A guard: without this, both tests above would pass on an empty matcher."""
    broken = ast.parse(
        "import marimo\n"
        "app = marimo.App()\n"
        "@app.cell\n"
        "def _():\n"
        "    import numpy as np\n"
        "    return (np,)\n"
        "@app.cell\n"
        "def _(plt):\n"
        "    return (np,)\n"
    )
    cells = _cells(broken)
    assert len(cells) == 2
    provided = set().union(*(_provided(c) for c in cells))
    consumed = set().union(*(_consumed(c) for c in cells))
    assert consumed - provided - BUILTINS == {"plt"}          # the aging.py failure
    assert _provided(cells[0]) & _provided(cells[1]) == {"np"}  # the model_explorer.py failure


def test_a_pass_through_is_not_a_second_definition() -> None:
    """The other half of the guard: re-returning a parameter must not count as a clash."""
    ok = ast.parse(
        "import marimo\n"
        "app = marimo.App()\n"
        "@app.cell\n"
        "def _():\n"
        "    import marimo as mo\n"
        "    return (mo,)\n"
        "@app.cell\n"
        "def _(mo):\n"
        "    import polars as pl\n"
        "    return mo, pl\n"
    )
    cells = _cells(ok)
    assert _provided(cells[0]) == {"mo"}
    assert _provided(cells[1]) == {"pl"}


def test_a_nested_functions_return_is_not_the_cells_return() -> None:
    """``emerson_cmv_hla.py`` defines a helper whose own ``return cmv, a02`` is not the cell's."""
    nested = ast.parse(
        "import marimo\n"
        "app = marimo.App()\n"
        "@app.cell\n"
        "def _():\n"
        "    def phenotypes(meta):\n"
        "        cmv = meta\n"
        "        return (cmv,)\n"
        "    return (phenotypes,)\n"
    )
    assert _provided(_cells(nested)[0]) == {"phenotypes"}
