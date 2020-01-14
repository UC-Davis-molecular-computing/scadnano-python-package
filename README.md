# scadnano-python-package

The scadnano Python module is a library for describing synthetic DNA nanostructures (e.g., DNA origami).

This module is used to write Python scripts outputting *.dna files readable by [scadnano](https://web.cs.ucdavis.edu/~doty/scadnano/), a web application useful for displaying and manually editing these structures. The purpose of this module is to help automate some of the task of creating DNA designs, as well as making large-scale changes to them that are easier to describe programmatically than to do by hand in scadnano.

## Installation

The scadnano Python package requires [Python version 3.7](https://www.python.org/downloads/) or later. 

There are two ways you can install it:


### pip 

1. If your Python installation does not already have pip installed, you may have to install it. 
Executing [this Python script](https://bootstrap.pypa.io/get-pip.py) should work; 
see also 
https://docs.python.org/3/installing/index.html 
or 
https://www.liquidweb.com/kb/install-pip-windows/.

2. Clone this repository by downloading [Git](https://git-scm.com/) and executing the following at the command line: 
    ```console
    git clone https://github.com/UC-Davis-molecular-computing/scadnano-python-package.git
    ```


3. Use pip to install the package from the local repository by executing the following at the command line from the same directory where you executed the git command above: 
    ```console
    pip install -e scadnano-python-package/
    ```

### download
As a simple alternative, you can download and place the following file(s) (located in the [scadnano/ subfolder](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/tree/master/scadnano)) in your PYTHONPATH (e.g., in the same directory as the scripts you are running):

* *required*: [scadnano.py](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/blob/master/scadnano/scadnano.py) 
* *optional*: [origami_rectangle.py](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/blob/master/scadnano/origami_rectangle.py); This can help create origami rectangles, but it is not necessary to use scadnano.

## Examples

Several example scripts are located in the [examples/](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/tree/master/examples) subfolder. Their output is contained in the [examples/output_designs/](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/tree/master/examples/output_designs) subfolder.

## Documentation

Online documentation of the package API is located here:
https://web.cs.ucdavis.edu/~doty/scadnano/docs/