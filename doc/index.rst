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
.. automodule:: scadnano.scadnano
   :members:

origami_rectangle 
=====================
.. automodule:: scadnano.origami_rectangle
   :members:

Interoperability - cadnano v2
=============================

Scadnano provides function to convert design to and from cadnano v2:

* :any:`DNADesign.from_cadnano_v2` will create a scadnano DNADesign from a ``cadnanov2`` json file. 
* :any:`DNADesign.export_cadnano_v2` will produce a ``cadnanov2`` json file from a scadnano design.

**Important**

All ``cadnanov2`` designs can be imported to scadnano. However **not all scadnano designs can be imported 
to cadnanov2**, to be importable to ``cadnanov2`` a scadnano design need to comply with the following points:

* The design cannot feature any :any:`Loopout` as it is not a concept that exists in ``cadnanov2``.
* Following ``cadnanov2`` conventions, helices with **even** number must have their scaffold going **forward** and helices with **odd** number **backward**.

Also note that maximum helices offsets can be altered in a ``scadnano`` to ``cadnanov2`` conversion as ``cadnanov2`` needs max offsets to be a multiple of 21 in the hex grid and 32 in the rectangular grid.
The conversion algorithm will choose the lowest multiple of 21 or 32 which fits the entire design.

The ``cadnanov2`` json format does not embed sequences hence they will be lost after conversion.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |br| raw:: html

   <br />

