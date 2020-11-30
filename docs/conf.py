

# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
import warnings
from pathlib import Path
from datetime import datetime
import sphinx_rtd_theme
#sys.path.insert(0, os.path.abspath('../py/scherlock'))
sys.path.insert(0, os.path.abspath('../py'))


HERE = Path(__file__).parent

# -- Project information -----------------------------------------------------

project = 'SCherlock'
copyright = '2020, Anton Björk, Johan Henriksson, Henlab'
author = 'Anton Björk, Johan Henriksson, Henlab'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.doctest',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx_autodoc_typehints',
    'sphinx_rtd_theme'
]



# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = 'bysource'
# autodoc_default_flags = ['members']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [('Params', 'Parameters')]
todo_include_todos = False
api_dir = HERE / 'api'  # function_images


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



intersphinx_mapping = dict(
    anndata=('https://anndata.readthedocs.io/en/stable/', None),
    bbknn=('https://bbknn.readthedocs.io/en/latest/', None),
    cycler=('https://matplotlib.org/cycler/', None),
    h5py=('http://docs.h5py.org/en/stable/', None),
    ipython=('https://ipython.readthedocs.io/en/stable/', None),
    leidenalg=('https://leidenalg.readthedocs.io/en/latest/', None),
    louvain=('https://louvain-igraph.readthedocs.io/en/latest/', None),
    matplotlib=('https://matplotlib.org/', None),
    networkx=('https://networkx.github.io/documentation/networkx-1.10/', None),
    numpy=('https://docs.scipy.org/doc/numpy/', None),
    pandas=('https://pandas.pydata.org/pandas-docs/stable/', None),
    python=('https://docs.python.org/3', None),
    scipy=('https://docs.scipy.org/doc/scipy/reference/', None),
    scvelo=('https://scvelo.readthedocs.io/', None),
    seaborn=('https://seaborn.pydata.org/', None),
    sklearn=('https://scikit-learn.org/stable/', None),
)

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
#html_theme_options = dict(navigation_depth=4, logo_only=True)  # Only show the logo
html_context = dict(
    display_github=True,  # Integrate GitHub
    github_user='henriksson-lab',  # Username
    github_repo='isocounter',  # Repo name
    github_version='master',  # Version
    conf_py_path='/docs/',  # Path in the checkout to the docs root
)
html_static_path = ['_static']
html_show_sphinx = False
#html_logo = '_static/img/Scanpy_Logo_BrightFG.svg'


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
