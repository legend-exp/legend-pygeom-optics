# Configuration file for the Sphinx documentation builder.
from __future__ import annotations

import sys
from pathlib import Path

from pkg_resources import get_distribution

sys.path.insert(0, Path(__file__).parents[2].resolve().as_posix())
# Add local extension directory.
sys.path.insert(0, (Path(__file__).parents[0] / "exts").as_posix())

project = "legend-pygeom-optics"
copyright = "The LEGEND Collaboration"
version = get_distribution("legend-pygeom-optics").version

extensions = [
    "sphinx.ext.githubpages",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx_copybutton",
    "myst_parser",
    "optics_plot_extension",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"
language = "python"

# Furo theme
html_theme = "furo"
html_theme_options = {
    "source_repository": "https://github.com/legend-exp/legend-pygeom-optics",
    "source_branch": "main",
    "source_directory": "docs/source",
}
html_title = f"{project} {version}"

autodoc_default_options = {
    "ignore-module-all": True,
    # ignore some common members from NamedTuples.
    "exclude-members": "_asdict, _fields, _field_defaults, _make, _replace",
}

# sphinx-napoleon
# enforce consistent usage of NumPy-style docstrings
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_ivar = True

# intersphinx
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
    "pint": ("https://pint.readthedocs.io/en/stable", None),
    "pyg4ometry": ("https://pyg4ometry.readthedocs.io/en/stable", None),
}  # add new intersphinx mappings here

# sphinx-autodoc
# Include __init__() docstring in class docstring
autoclass_content = "both"
autodoc_typehints = "description"
autodoc_typehints_description_target = "all"
autodoc_typehints_format = "short"
# mock pyg4ometry that will sometimes lead to build failures.
autodoc_mock_imports = ["pyg4ometry"]
