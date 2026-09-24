"""Sphinx configuration for the vdjtools documentation."""

# Read the version from the installed package rather than repeating it here: this file said
# 3.6.1 through six minor releases, so every docs page carried a wrong version in its header.
from vdjtools import __version__

project = "vdjtools"
author = "ISALGO laboratory"
copyright = "2026, ISALGO laboratory"
version = release = __version__

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "nbsphinx",
]

# The compiled _core ext is installed in the docs build env, so vdjtools imports; the
# heavy parent deps (only imported by the model / overlap subpackages) are mocked.
autodoc_mock_imports = ["arda", "vdjmatch", "seqtree"]
autodoc_typehints = "description"
autodoc_member_order = "bysource"

# Render napoleon ``Attributes:`` sections as :ivar: fields (not standalone py:attribute
# objects), so a dataclass's Attributes docstring does not duplicate its autodoc'd fields.
napoleon_use_ivar = True

intersphinx_mapping = {"python": ("https://docs.python.org/3", None)}

templates_path = ["_templates"]
exclude_patterns = ["_build", "**.ipynb_checkpoints"]
html_static_path = ["_static"]
html_css_files = ["custom.css"]

html_theme = "pydata_sphinx_theme"
html_title = f"vdjtools {release}"
html_theme_options = {
    # Version shown in the navbar brand on every page (no image logo -> text brand).
    "logo": {"text": f"vdjtools {release}"},
    "github_url": "https://github.com/antigenomics/vdjtools",
    "navigation_with_keys": True,
    "show_prev_next": False,
    # **The sidebar was a flat wall of page titles** in source order, with no grouping and nothing
    # marking where the reader is -- so it was scenery, not navigation, and `signature` sat in it
    # indistinguishable from everything else. `index.rst` now carries captioned toctrees; these
    # options make the theme render that structure rather than flatten it. Same settings as
    # `seqtree` and `mhcmatch`.
    "show_nav_level": 2,        # open each caption's pages, do not collapse to the caption alone
    "navigation_depth": 3,      # let a page's own sections show under it
    "collapse_navigation": False,
    "header_links_before_dropdown": 4,
    "show_toc_level": 2,        # right-hand "On this page": subsections too, not just top level
}

# The stock sidebar renders only children of the current top-level page, and these pages are all
# top-level siblings -- so it renders nothing. `site-nav` renders the full captioned tree.
html_sidebars = {
    "**": ["site-nav"],
    "index": [],
}
nbsphinx_execute = "never"
