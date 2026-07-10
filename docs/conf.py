"""Sphinx configuration for the vdjtools documentation."""

project = "vdjtools"
author = "ISALGO laboratory"
copyright = "2026, ISALGO laboratory"
version = release = "2.0.0"

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

intersphinx_mapping = {"python": ("https://docs.python.org/3", None)}

templates_path = ["_templates"]
exclude_patterns = ["_build", "**.ipynb_checkpoints"]

html_theme = "pydata_sphinx_theme"
html_title = "vdjtools"
html_theme_options = {
    "github_url": "https://github.com/antigenomics/vdjtools",
    "navigation_with_keys": True,
}
nbsphinx_execute = "never"
