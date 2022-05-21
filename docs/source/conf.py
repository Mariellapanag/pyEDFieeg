# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'pyEDFieeg'
copyright = '2022, Mariella Panagiotopoulou'
author = 'Mariella Panagiotopoulou'

# The full version, including alpha/beta/rc tags
release = '0.0.1'


# General
extensions = [
    "myst_parser",
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx_autodoc_typehints",
    "sphinx.ext.viewcode",
    "sphinx.ext.doctest",
]
exclude_patterns = ["_build"]

# HTML output
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "navigation_depth": 2,
}
html_context = {
    "conf_py_path": "/docs",
    "source_suffix": ".md",
}
html_title = "torchtime"
html_static_path = ["_static"]
pygments_style = "friendly"
html_show_sphinx = False
