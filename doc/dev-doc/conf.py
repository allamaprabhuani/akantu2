# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config
__copyright__ = (
    "Copyright (©) 2020-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import glob
import shutil
import subprocess
import sys
import jinja2


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Number figures
numfig = True

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinxcontrib.bibtex",
    "breathe",
    "myst_parser",
#    "sphinx_gallery.gen_gallery",
    "sphinx_copybutton",
]

read_the_docs_build = os.environ.get("READTHEDOCS", None) == "True"
cmake_configure = os.environ.get("RUNNING_IN_CMAKE", None) == "True"
if read_the_docs_build:
    akantu_path = "."
    akantu_source_path = "../../"
elif cmake_configure:
    akantu_path = "@CMAKE_CURRENT_BINARY_DIR@"
    akantu_source_path = "@CMAKE_SOURCE_DIR@"
else:  # most probably running by hand
    akantu_path = "build-doc"
    akantu_source_path = "../../"
    try:
        os.mkdir(akantu_path)
    except FileExistsError:
        pass


if akantu_path == "@" + "CMAKE_CURRENT_BINARY_DIR" + "@":  # Concatenation is to avoid cmake to replace it
    raise Exception("Something went really wrong")

sys.path.insert(0, akantu_source_path)

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "CMakeLists.txt",
    "manual/appendix/elements.rst",
    "manual/appendix/material-parameters.rst",
    "manual/constitutive-laws.rst",
    "manual/new-constitutive-laws.rst",
    "examples/README.rst",
    "examples/c++/README.rst",
    "examples/c++/contact_mechanics_model/README.rst",
    "examples/c++/solid_mechanics_model/README.rst",
    "examples/c++/solid_mechanics_model/boundary_conditions/README.rst",
    "examples/c++/solid_mechanics_model/static/README.rst",
    "examples/c++/solid_mechanics_model/explicit/README.rst",
    "examples/c++/solid_mechanics_cohesive_model/README.rst",
    "examples/python/README.rst",
]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

primary_domain = "cpp"
highlight_language = "cpp"

bibtex_bibfiles = ["manual/manual-bibliography.bib"]

# -- Project information -----------------------------------------------------

project = "Akantu"
copyright = (
    "2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    + " Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
author = "Nicolas Richart"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
if read_the_docs_build:
    html_theme = "default"
else:
    html_theme = "sphinx_rtd_theme"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_logo = "_static/logo_only_akantu.svg"
# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}
html_sidebars = {
    "**": [
        "relations.html",  # needs 'show_related': True theme option to display
        "searchbox.html",
    ]
}

math_eqref_format = "Eq. {number}"

# MathJax configuration
mathjax3_config = {
    "tex": {
        "macros": {
            "st": ["\\mathrm{#1}", 1],
            "mat": ["\\mathbf{#1}", 1],
            "half": "\\frac{1}{2}",
        },
        "packages": {"[+]": ["ams"]},
    },
    "loader": {"load": ["[tex]/ams"]},
}

# -- Options for HTMLHelp output ---------------------------------------------
# Output file base name for HTML help builder.
htmlhelp_basename = "Akantudoc"


# -- Options for LaTeX output ------------------------------------------------
latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    #
    "preamble": r"""\usepackage{amsmath}""",
    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, "Akantu.tex", "Akantu Documentation",
     "Nicolas Richart", "manual"),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "akantu", "Akantu Documentation", [author], 1)]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "Akantu",
        "Akantu Documentation",
        author,
        "Akantu",
        "One line description of project.",
        "Miscellaneous",
    ),
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
epub_identifier = ""

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ["search.html"]


# -- Extension configuration -------------------------------------------------
j2_args = {}

if read_the_docs_build or not cmake_configure:
    j2_template_path = "."
elif cmake_configure:
    j2_template_path = "@CMAKE_CURRENT_SOURCE_DIR@"
    os.makedirs(os.path.join(akantu_path, "_static"), exist_ok=True)
    shutil.copyfile(
        os.path.join("@CMAKE_CURRENT_SOURCE_DIR@", html_logo),
        os.path.join(akantu_path, html_logo),
    )

j2_args = {
    "akantu_source_path": akantu_source_path,
}

print(akantu_path)
j2_env = jinja2.Environment(
    loader=jinja2.FileSystemLoader(j2_template_path),
    undefined=jinja2.DebugUndefined
)

j2_template = j2_env.get_template("akantu.dox.j2")

with open(os.path.join(akantu_path, "akantu.dox"), "w") as fh:
    fh.write(j2_template.render(j2_args))

subprocess.run(["doxygen", "-q", "akantu.dox"], cwd=akantu_path)

# print("akantu_path = '{}'".format(akantu_path))
breathe_projects = {"Akantu": os.path.join(akantu_path, "xml")}
breathe_default_project = "Akantu"
breathe_default_members = ("members", "undoc-members")
breathe_implementation_filename_extensions = [".c", ".cc", ".cpp", ".hh"]
breathe_show_enumvalue_initializer = True
breathe_debug_trace_directives = False
breathe_short_warning = True

# -- Gallery ------------------------------------------------------------------
# sphinx_gallery_conf = {
#     'examples_dirs': os.path.join(akantu_source_path, 'examples'),
#     'gallery_dirs': os.path.join(akantu_source_path, 'doc', 'dev-doc', 'auto_examples'),
#     'download_all_examples': False,
#     'plot_gallery': 'False',
#     'only_warn_on_example_error': True,
#     'log_level': {'backreference_missing': 'debug'},
# }


# -- Options for intersphinx extension ----------------------------------------
intersphinx_mapping = {
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference", None),
}
