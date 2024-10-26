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

Scadnano provides function to convert design to and from cadnano v2:

* :meth:`scadnano.DNADesign.from_cadnano_v2` will create a scadnano DNADesign from a ``cadnanov2`` json file.
* :meth:`scadnano.DNADesign.export_cadnano_v2` will produce a ``cadnanov2`` json file from a scadnano design.

**Important**

All ``cadnanov2`` designs can be imported to scadnano. However **not all scadnano designs can be imported 
to cadnanov2**; to be exportable to ``cadnanov2`` a scadnano design need to comply with the following constraints:

* The design cannot feature any :class:`Loopout` or :class:`Extension`, since these are not concepts that exist in
  ``cadnanov2``.
* Following ``cadnanov2`` conventions, helices with **even** number must have their scaffold going **forward** and
  helices with **odd** number **backward**.
* If you use paranemic crossovers (i.e. crossovers where the domains before and after the crossover are in the same
  direction), the helices' row number (i.e. not the helices' indexes but their y coordinate) of the domains must have
  the same parity, meaning both even or both odd, for example rows 0 and 2, or 1 and 3, but not 0 and 1.

Also note that maximum helices offsets can be altered in a ``scadnano`` to ``cadnanov2`` conversion as ``cadnanov2``
needs max offsets to be a multiple of 21 in the hex grid and 32 in the rectangular grid. The conversion algorithm will
choose the lowest multiple of 21 or 32 which fits the entire design.

The ``cadnanov2`` json format does not embed sequences hence they will be lost after conversion.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |br| raw:: html

   <br />

