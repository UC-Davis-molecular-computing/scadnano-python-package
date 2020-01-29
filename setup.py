#!/usr/bin/env python

from distutils.core import setup

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='scadnano',
      packages=['scadnano'],
      version='0.1.1',
      license='MIT',
      description="Python scripting library for generating designs readable by scadnano.",
      author="David Doty",
      author_email="doty@ucdavis.edu",
      url="https://github.com/UC-Davis-molecular-computing/scadnano-python-package",
      download_url = 'https://github.com/UC-Davis-molecular-computing/scadnano-python-package/archive/v0.1.0.zip', 
      long_description = long_description,
      long_description_content_type = 'text/markdown',
    )
