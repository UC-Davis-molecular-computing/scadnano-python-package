#!/usr/bin/env python

# from distutils.core import setup

# got some ideas from here: https://medium.com/@joel.barmettler/how-to-upload-your-python-package-to-pypi-65edc5fe9c56
#
# But apparently setuptools is the replacement for distutils, and distutils was causing problems such as
# not including the README.md file and not formatting it as Markdown on PyPI
# https://setuptools.readthedocs.io/en/latest/setuptools.html

from setuptools import setup


# import scadnano.scadnano_version as sv


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


__version__ = extract_version('scadnano/scadnano.py')

# read the contents of your README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='scadnano',
      packages=['scadnano'],
      version=__version__,
      # version='0.8.0',
      download_url=f'https://github.com/UC-Davis-molecular-computing/scadnano-python-package/archive/v{__version__}.zip',
      # download_url=f'https://github.com/UC-Davis-molecular-computing/scadnano-python-package/archive/v0.7.0.zip',
      license='MIT',
      description="Python scripting library for generating designs readable by scadnano.",
      author="David Doty",
      author_email="doty@ucdavis.edu",
      url="https://github.com/UC-Davis-molecular-computing/scadnano-python-package",
      long_description=long_description,
      long_description_content_type='text/markdown; variant=GFM',
      python_requires='>=3.6',
      install_requires=[
        'xlwt',
        'dataclasses>=0.6; python_version < "3.7"',
        'tabulate',
      ],
      tests_require=['xlrd'],
)
