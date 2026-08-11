"""Sphinx configuration for the CpGtools documentation."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version as package_version


# -- Project information -----------------------------------------------------

project = "CpGtools"
author = "Liguo Wang"
copyright = "2024-2026, Liguo Wang"

try:
    release = package_version("cpgtools")
except PackageNotFoundError:
    # Fallback used when CpGtools is not installed in the documentation
    # environment. Update this value when preparing a release if needed.
    release = "3.0.0"

version = ".".join(release.split(".")[:2])


# -- General configuration ---------------------------------------------------

needs_sphinx = "7.2"

extensions = [
    "sphinx_rtd_theme",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

source_suffix = {
    ".rst": "restructuredtext",
}

root_doc = "index"
language = "en"
pygments_style = "sphinx"


# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_title = "CpGtools Documentation"
html_short_title = "CpGtools"
html_logo = "_static/CpGtools_logo3.png"
html_static_path = ["_static"]
html_last_updated_fmt = "%B %d, %Y"
html_show_sphinx = False
html_show_copyright = True

html_theme_options = {
    "logo_only": True,
    "display_version": True,
    "prev_next_buttons_location": "bottom",
    "style_external_links": False,
    "collapse_navigation": True,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "includehidden": True,
    "titles_only": False,
}


# -- Options for HTML Help output -------------------------------------------

htmlhelp_basename = "CpGtoolsdoc"


# -- Options for LaTeX output ------------------------------------------------

latex_elements: dict[str, str] = {}

latex_documents = [
    (
        "index",
        "CpGtools.tex",
        "CpGtools Documentation",
        author,
        "manual",
    ),
]


# -- Options for manual page output -----------------------------------------

man_pages = [
    (
        "index",
        "cpgtools",
        "CpGtools Documentation",
        [author],
        1,
    ),
]


# -- Options for Texinfo output ---------------------------------------------

texinfo_documents = [
    (
        "index",
        "CpGtools",
        "CpGtools Documentation",
        author,
        "CpGtools",
        "Python package for DNA methylation data analysis.",
        "Bioinformatics",
    ),
]


def setup(app) -> None:
    """Register project-specific static assets."""
    app.add_css_file("custom.css")
