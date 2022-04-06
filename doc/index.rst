.. scadnano documentation master file, created by
   sphinx-quickstart on Tue Jul 16 10:14:04 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

scadnano documentation
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

scadnano 
=====================
.. automodule:: scadnano
   :members:

origami_rectangle 
=====================
.. automodule:: origami_rectangle
   :members:

Interoperability - cadnano v2
=============================

The scadnano Python package provides functionality to convert design to and from cadnano v2:

* :meth:`scadnano.Design.from_cadnano_v2` will create a scadnano DNADesign from a ``cadnanov2`` json file.
* :meth:`scadnano.Design.to_cadnano_v2_json` will produce a ``cadnanov2`` json format from a scadnano design 
  (returning the result as a string).
* :meth:`scadnano.Design.write_cadnano_v2_file` will produce a ``cadnanov2`` json file from a scadnano design 
  and write it to disk.

Import and export are also supported in the web interface (https://scadnano.org) 
via the menus "File→Import cadnano v2" and "Export→cadnano v2".

**Important**

All ``cadnanov2`` designs can be imported to scadnano. However **not all scadnano designs can be exported 
to cadnanov2**, to be exportable to ``cadnanov2`` a scadnano design need to comply with the following points:

* The design cannot feature any :class:`Loopout` as it is not a concept that exists in ``cadnanov2``.
* Following ``cadnanov2`` conventions, helices with **even** number must have their scaffold going **forward** 
  and helices with **odd** number **backward**.

Maximum helix offsets may be altered in a ``scadnano`` to ``cadnanov2`` conversion, because 
``cadnanov2`` needs max offsets to be a multiple of 21 in the honeycomb grid and a multiple of 32 in the square grid.
The conversion algorithm will choose the lowest multiple of 21 or 32 that fits the entire design.
The hex grid and none grid are not supported.

The ``cadnanov2`` json format does not embed DNA sequences, so they will be lost in the cadnano export.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |br| raw:: html

   <br />

