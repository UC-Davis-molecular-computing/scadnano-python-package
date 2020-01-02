# scadnano-python-package

## Installation

The scadnano Python package requires [Python version 3.7](https://www.python.org/downloads/) or later. 


### pip 

1. If your Python installation does not already have pip installed, you may have to install it: see https://docs.python.org/3/installing/index.html or https://www.liquidweb.com/kb/install-pip-windows/

2. Clone this repository by downloading [git](https://git-scm.com/) and executing the following at the command line: `git clone https://github.com/UC-Davis-molecular-computing/scadnano-python-package.git`

3. Use pip to install the package from the local repository by executing the following at the command line from the same directory where you executed the git command above: `pip install -e scadnano-python-package/`

### copy files
As a simple alternative, you can download and place the following files (from the scadnano/ subdirectory of the [repository](https://github.com/UC-Davis-molecular-computing/scadnano-python-package)) in your PYTHONPATH (e.g., in the same directory as the scripts you are running):

* [scadnano.py](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/blob/master/scadnano/scadnano.py)
* [m13.py](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/blob/master/scadnano/m13.py)
* [json_utils.py](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/blob/master/scadnano/json_utils.py)
* [origami_rectangle.py](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/blob/master/scadnano/origami_rectangle.py)

## Examples

Several example scripts are located in the examples/ subfolder. Their output is contained in the examples/output_designs/ subfolder.

## Documentation

Online documentation of the package API is located here:
https://web.cs.ucdavis.edu/~doty/scadnano/docs/