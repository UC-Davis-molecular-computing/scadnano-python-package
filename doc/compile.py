""" Take index.rst.in as input and populate with classes of scadnano.py
"""

import scadnano as sc 

for elem in dir(sc):
    if 