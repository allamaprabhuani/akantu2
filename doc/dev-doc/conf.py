# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import shutil
import jinja2
import git
import re
import subprocess

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
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinxcontrib.bibtex',
    'breathe',
]

read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'
if read_the_docs_build:
    akantu_path = "."
    akantu_source_path = "../../"
else:
    akantu_path = "@CMAKE_CURRENT_BINARY_DIR@"
    akantu_source_path = "@CMAKE_SOURCE_DIR@"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['CMakeLists.txt',
                    'manual/appendix/elements.rst',
                    'manual/appendix/material-parameters.rst',
                    'manual/constitutive-laws.rst',
                    'manual/new-constitutive-laws.rst']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

primary_domain = 'cpp'
highlight_language = 'cpp'

bibtex_bibfiles = ['manual/manual-bibliography.bib']

# -- Project information -----------------------------------------------------

project = 'Akantu'
copyright = '2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)' + \
    ' Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)'
author = 'Nicolas Richart'

with open(os.path.join(akantu_source_path, 'VERSION'), 'r') as fh:
    version_file = fh.readlines()
    file_release = version_file[0].strip()

try:
    tag_prefix = 'v'
    git_repo = git.Repo(akantu_source_path)

    git_describe = git_repo.git.describe('--tags', '--dirty', '--always',
                                         '--long',
                                         '--match', '{}*'.format(tag_prefix))

    print("GIT Describe: {}".format(git_describe))

    # git describe to PEP404 version
    describe_matches = re.search(
        (r'^{}(?P<version>.+?)' +
         r'(?:-(?P<distance>\d+)-g(?P<sha>[0-9a-f]+)' +
         r'(?:-(?P<dirty>dirty))?)?$').format(tag_prefix),
        git_describe)

    if describe_matches:
        describe_matches = describe_matches.groupdict()

        release = describe_matches['version']
        if describe_matches['distance']:
            release += '.' if '+' in release else '+'
            release += '{distance}.{sha}'.format(**describe_matches)
            if describe_matches['dirty']:
                release += '.dirty'
    else:
        count = git_repo.git.rev_list('HEAD', '--count')
        describe_matches = re.search(
            (r'^(?P<sha>[0-9a-f]+)' +
             r'(?:-(?P<dirty>dirty))?$').format(tag_prefix),
            git_describe).groupdict()
        release = '{}.{}+{}'.format(file_release, count,
                                    describe_matches['sha'])

except git.InvalidGitRepositoryError:
    release = file_release

version = re.sub(r'^([0-9]+)\.([0-9+]).*',
                 r'\1.\2',
                 release)

print('Release: {} - Version: {}'.format(release, version))

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
if read_the_docs_build:
    html_theme = 'default'
else:
    html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_logo = '_static/logo_only_akantu.svg'
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
    '**': [
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
    ]}

math_eqref_format = "Eq. {number}"

# MathJax configuration
if not read_the_docs_build:
    mathjax_config = {
        'extensions': [
            "tex2jax.js",
            "siunitx.js"
        ],
        'TeX': {
            'Macros': {
                'st': [r'\mathrm{#1}', 1],
                'mat': [r'\mathbf{#1}', 1],
                'half': [r'\frac{1}{2}', 0],
            },
            'extensions': ["AMSmath.js", "AMSsymbols.js", "sinuitx.js"],
        },
    }
else:
    mathjax3_config = {
        'tex': {
            'macros': {
                'st': [r'\mathrm{#1}', 1],
                'mat': [r'\mathbf{#1}', 1],
                'half': [r'\frac{1}{2}', 0],
            },
            'packages': ['base', 'ams'],
        },
        'loader': {
            'load': ['[tex]/ams']
        },
    }


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'Akantudoc'


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
    'preamble': r'''\usepackage{amsmath}''',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'Akantu.tex', 'Akantu Documentation',
     'Nicolas Richart', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'akantu', 'Akantu Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'Akantu', 'Akantu Documentation',
     author, 'Akantu', 'One line description of project.',
     'Miscellaneous'),
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
epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']


# -- Extension configuration -------------------------------------------------
j2_args = {}

if read_the_docs_build:
    j2_template_path = '.'
else:
    j2_template_path = '@CMAKE_CURRENT_SOURCE_DIR@'
    os.makedirs(os.path.join(akantu_path, '_static'), exist_ok=True)
    shutil.copyfile(
        os.path.join('@CMAKE_CURRENT_SOURCE_DIR@', html_logo),
        os.path.join(akantu_path, html_logo))

j2_args = {
    'akantu_source_path': akantu_source_path,
    'akantu_version': version.replace('v', ''),
}

print(akantu_path)
j2_env = jinja2.Environment(
    loader=jinja2.FileSystemLoader(j2_template_path),
    undefined=jinja2.DebugUndefined)

j2_template = j2_env.get_template('akantu.dox.j2')

with open(os.path.join(akantu_path, 'akantu.dox'), 'w') as fh:
    fh.write(j2_template.render(j2_args))

subprocess.run(['doxygen', 'akantu.dox'],
               cwd=akantu_path)

# print("akantu_path = '{}'".format(akantu_path))
breathe_projects = {"Akantu": os.path.join(akantu_path, "xml")}
breathe_default_project = "Akantu"
breathe_default_members = ('members', 'undoc-members')
breathe_implementation_filename_extensions = ['.c', '.cc', '.cpp']
breathe_show_enumvalue_initializer = True

# -- Options for intersphinx extension ---------------------------------------

intersphinx_mapping = {
    'numpy': ('https://docs.scipy.org/doc/numpy/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference', None),
}
