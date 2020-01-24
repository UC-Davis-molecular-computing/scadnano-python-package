# scadnano-python-package

The scadnano Python module is a library for describing synthetic DNA nanostructures (e.g., DNA origami).

This module is used to write Python scripts outputting *.dna files readable by [scadnano](https://web.cs.ucdavis.edu/~doty/scadnano/), a web application useful for displaying and manually editing these structures. The purpose of this module is to help automate some of the task of creating DNA designs, as well as making large-scale changes to them that are easier to describe programmatically than to do by hand in scadnano.

## Installation

The scadnano Python package requires [Python version 3.7](https://www.python.org/downloads/) or later. 

There are two ways you can install the scadnano Python package:


1. pip 

    Use [pip](https://pypi.org/project/pip/) to install the package by executing the following at the command line:
    ```console
    pip install scadnano
    ```

    If your Python installation does not already have pip installed, you may have to install it. 
    Executing [this Python script](https://bootstrap.pypa.io/get-pip.py) should work; 
    see also 
    https://docs.python.org/3/installing/index.html 
    or 
    https://www.liquidweb.com/kb/install-pip-windows/.

2. download

    As a simple alternative, you can download and place the following file(s) (located in the [scadnano/](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/tree/master/scadnano) subfolder)
    in your PYTHONPATH (e.g., in the same directory as the scripts you are running):

    * *required*: [scadnano.py](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/blob/master/scadnano/scadnano.py) 
    * *optional*: [origami_rectangle.py](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/blob/master/scadnano/origami_rectangle.py); This can help create origami rectangles, but it is not necessary to use scadnano.

## Documentation

Online documentation of the package API is located here:
https://web.cs.ucdavis.edu/~doty/scadnano/docs/

## Tutorial

A [tutorial](tutorial.md) is available.

## Examples

Several example scripts are located in the 
[examples/](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/tree/master/examples) subfolder. 
Their output is contained in the 
[examples/output_designs/](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/tree/master/examples/output_designs) subfolder.
