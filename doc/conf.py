# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
import tomllib

project = 'TEDS'
copyright = 'SRON'
author = 'TANGO team'

pyproject_data = tomllib.load(open('../pyproject.toml', 'rb'))
version = pyproject_data['projects']['version']

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.napoleon']

templates_path = ['_templates']
exclude_patterns = ['_build']

latex_elements = {
  'extraclassoptions': 'openany'
}

# Napoleon settings
sys.path.insert(0, os.path.abspath('..'))
napoleon_numpy_docstring = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'piccolo_theme'
html_static_path = ['_static']
html_title = 'TEDS'
html_short_title = 'TEDS'
