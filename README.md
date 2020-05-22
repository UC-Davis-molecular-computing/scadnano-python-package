# scadnano-python-package

![Python package](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/workflows/Python%20package/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/scadnano-python-package/badge/?version=latest)](https://scadnano-python-package.readthedocs.io/en/latest/?badge=latest)


The scadnano Python module is a library for describing synthetic DNA nanostructures (e.g., DNA origami).

This module is used to write Python scripts outputting `.dna` files readable by [scadnano](https://scadnano.org), a web application useful for displaying and manually editing these structures. The purpose of this module is to help automate some of the task of creating DNA designs, as well as making large-scale changes to them that are easier to describe programmatically than to do by hand in scadnano.

## Reporting issues

Please report issues in the web interface at the [scadnano web interface GitHub repository](https://github.com/UC-Davis-molecular-computing/scadnano/issues), and report issues in the Python scripting library at the [scadnano Python package GitHub repository](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues).

## Installation

The scadnano Python package requires Python version 3.7 or later. If you do not have that version (or later) of Python installed, follow [this link](https://www.python.org/downloads/) to install it.

Once Python is installed, there are two ways you can install the scadnano Python package:


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

    * *required*: [scadnano.py](https://raw.githubusercontent.com/UC-Davis-molecular-computing/scadnano-python-package/master/scadnano/scadnano.py) 
    * *optional*: [modifications.py](https://raw.githubusercontent.com/UC-Davis-molecular-computing/scadnano-python-package/master/scadnano/modifications.py); This contains some common DNA modifications such as biotin and Cy3. 
    * *optional*: [origami_rectangle.py](https://raw.githubusercontent.com/UC-Davis-molecular-computing/scadnano-python-package/master/scadnano/origami_rectangle.py); This can help create origami rectangles, but it is not necessary to use scadnano.

    Unfortunately, the scadnano package depends on the Python package [xlwt](https://pypi.org/project/xlwt/), so xlwt must be installed. To install, type `pip install xlwt` at the command line. If this does not work, then you can comment out the line `import xlwt` in `scadnano.py`, but then you will not be able to call the method [`DNADesign.write_idt_plate_excel_file()`](https://scadnano-python-package.readthedocs.io/#scadnano.scadnano.DNADesign.write_idt_plate_excel_file) to export an Excel file with DNA sequences.





## Documentation

Online documentation of the package API is located here:
https://scadnano-python-package.readthedocs.io





## Tutorial

A [tutorial](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/blob/master/tutorial/tutorial.md) shows how to create a "standard" 24-helix DNA origami rectangle using the scripting library.





## Examples

Several example scripts are located in the 
[examples/](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/tree/master/examples) subfolder. 
Their output is contained in the 
[examples/output_designs/](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/tree/master/examples/output_designs) subfolder.
