# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

import sphinx.util.inspect as inspect
import sphinx.ext.autodoc as auto
import os
import sys



sys.path.insert(0, os.path.abspath('../scadnano'))
# print(sys.path)

# from scadnano.scadnano_version import current_version
# this is ugly, but appears to be standard practice:
# https://stackoverflow.com/questions/17583443/what-is-the-correct-way-to-share-package-version-with-setup-py-and-the-package/17626524#17626524
# __version__ = open("../scadnano/_version.py").readlines()[-1].split()[-1].strip("\"'")

# this is ugly, but appears to be standard practice:
# https://stackoverflow.com/questions/17583443/what-is-the-correct-way-to-share-package-version-with-setup-py-and-the-package/17626524#17626524
def extract_version(filename: str):
    with open(filename) as f:
        lines = f.readlines()
    version_comment = '# version line; WARNING: do not remove or change this line or comment'
    for line in lines:
        if version_comment in line:
            idx = line.index(version_comment)
            line_prefix = line[:idx]
            parts = line_prefix.split('=')
            stripped_parts = [part.strip() for part in parts]
            version_str = stripped_parts[-1].replace('"', '')
            return version_str
    raise AssertionError(f'could not find version in {filename}')

__version__ = extract_version('../scadnano/scadnano.py')


# Type "make html" at the command line to generate the documentation.


# -- Project information -----------------------------------------------------

project = 'scadnano'
copyright = '2019, David Doty'
author = 'David Doty'

# The full version, including alpha/beta/rc tags
# release = '0.7.0'
version = __version__
release = __version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = "sphinx_rtd_theme"
# html_theme = 'alabaster'
# html_theme = "classic"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

# use order in source rather than alphabetical order
autodoc_member_order = 'bysource'

# intersphinx_mapping = {'python': ('https://docs.python.org/3', None)}

# removes constant values from documentation that are longer than 100 characters.
# taken from
# https://stackoverflow.com/questions/25145817/ellipsis-truncation-on-module-attribute-value-in-sphinx-generated-documentatio/25163963#25163963


# from sphinx.ext.autodoc import DataDocumenter, ModuleLevelDocumenter, SUPPRESS
# from sphinx.util.inspect import safe_repr

length_limit = 50


def add_directive_header(self, sig):
    auto.ModuleLevelDocumenter.add_directive_header(self, sig)
    if not self.options.annotation:
        try:
            # objrepr = inspect.safe_repr(self.object)
            objrepr = inspect.object_description(self.object)

            # PATCH: truncate the value if longer than length_limit characters
            if len(objrepr) > length_limit:
                objrepr = objrepr[:length_limit] + "..."

        except ValueError:
            pass
        else:
            self.add_line(u'   :annotation: = ' + objrepr, '<autodoc>')
    elif self.options.annotation is auto.SUPPRESS:
        pass
    else:
        self.add_line(u'   :annotation: %s' % self.options.annotation,
                      '<autodoc>')


auto.DataDocumenter.add_directive_header = add_directive_header


# https://stackoverflow.com/questions/56336234/build-fail-sphinx-error-contents-rst-not-found
master_doc = 'index'
