project = "pylibecpint"
copyright = "2025, Robert A. Shaw"
author = "Robert A. Shaw"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
}

autodoc_member_order = "bysource"

html_theme = "furo"

exclude_patterns = ["_build"]
