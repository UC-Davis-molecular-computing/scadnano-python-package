# scadnano-python-package

![Python package](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/workflows/Python%20package/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/scadnano-python-package/badge/?version=latest)](https://scadnano-python-package.readthedocs.io/en/latest/?badge=latest)

The scadnano Python module is a library for describing synthetic DNA nanostructures (e.g., DNA origami).


If you find scadnano useful in a scientific project, please cite its associated paper:

> scadnano: A browser-based, scriptable tool for designing DNA nanostructures.  
  David Doty, Benjamin L Lee, and Tristan Stérin.  
  DNA 2020: *Proceedings of the 26th International Conference on DNA Computing and Molecular Programming*  
  [ [paper](https://arxiv.org/abs/2005.11841) | [BibTeX](https://web.cs.ucdavis.edu/~doty/papers/scadnano.bib) ]


## Overview

This module is used to write Python scripts outputting `.dna` files readable by [scadnano](https://scadnano.org), a web application useful for displaying and manually editing these structures. The purpose of this module is to help automate some of the task of creating DNA designs, as well as making large-scale changes to them that are easier to describe programmatically than to do by hand in scadnano.

Early versions of this project didn't have well-defined versions. However, we will try to announce breaking changes (and possibly new features) under the [GitHub releases page](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases). The version numbers in this Python library repo and the [web interface repo](https://github.com/UC-Davis-molecular-computing/scadnano/releases) won't always advance at the same time. 

Following [semantic versioning](https://semver.org/), version numbers are major.minor.patch, i.e., version 0.9.2 has minor version number 9. Prior to version 1.0.0, when a breaking change is made, this will increment the minor version (for example, going from 0.9.4 to 0.10.0). After version 1.0.0, breaking changes will increment the major version.





## Reporting issues

Please report issues in the web interface at the [scadnano web interface GitHub repository](https://github.com/UC-Davis-molecular-computing/scadnano/issues), and report issues in the Python scripting library at the [scadnano Python package GitHub repository](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues).






## Installation

The scadnano Python package requires Python version 3.7 or later. If you do not have that version (or later) of Python installed, follow [this link](https://www.python.org/downloads/) to install it. To check your current version of Python, open a command line and type

```
python --version
```

If you are using Python 3.6 and do not wish to upgrade, you can install a package providing the required features.
There is a Python 3.7 module used by scadnano, which is not present in 3.6: the
[dataclasses](https://docs.python.org/3/library/dataclasses.html) module.
If you have at least Python version 3.6, you can install a backported library providing this module: 
[dataclasses backport](https://pypi.org/project/dataclasses/)

Once Python is installed (and the dataclasses backport if you have Python version 3.6), there are two ways you can install the scadnano Python package:


1. pip (recommended)

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

    First, check your version of `pip` by typing 
    ```
    pip --version
    ```
    It should say something like 
    ```
    pip 19.3.1 from ...lib\site-packages\pip (python 3.8)
    ```
    If the version of Python at the end is Python 3.7 or higher, you are good. If it is version 2.7 or lower, type
    ```
    pip3 --version
    ```
    If that works and shows Python 3.7 or higher, you are good, but you should type `pip3` in the subsequent instructions instead of `pip`.
    
    If it shows Python 3.6, install the [dataclasses backport module](https://pypi.org/project/dataclasses/) via
    ```
    pip install dataclasses
    ```
    If it shows Python 3.5 or lower, then you will need to upgrade your Python version (recommended Python 3.7 or higher).

    
2. download

    *Note:* If you are reading this on the PyPI website, the links below won't work. They are relative links intended to be read on the [GitHub README page](https://github.com/UC-Davis-molecular-computing/scadnano-python-package#readme).

    As a simple alternative (in case you run into trouble using pip), you can download and place the following files (located in the [scadnano/](scadnano) subfolder) in your PYTHONPATH (e.g., in the same directory as the scripts you are running). To download them, right-click on "Raw" near the top and select (in Chrome or Firefox) "Save link as...". 

    * *required*: [scadnano.py](scadnano/scadnano.py) 
    * *optional*: [modifications.py](scadnano/modifications.py); This contains some common DNA modifications such as biotin and Cy3. 
    * *optional*: [origami_rectangle.py](scadnano/origami_rectangle.py); This can help create origami rectangles, but it is not necessary to use scadnano.
    * *optional*: [_version.py](scadnano/_version.py) This ensures that the current version number is written into any `.dna` files written by the library; otherwise it may be out of date. (Which should not matter for the most part.)
    
    The scadnano package uses the Python package [xlwt](https://pypi.org/project/xlwt/) to write Excel files, so xlwt must be installed in order to call the method [`DNADesign.write_idt_plate_excel_file()`](https://scadnano-python-package.readthedocs.io/#scadnano.DNADesign.write_idt_plate_excel_file) to export an Excel file with DNA sequences. To install xlwt, type `pip install xlwt` at the command line. (If you instead use pip to install the scadnano package, xlwt will be automatically installed.)







## Example

Consider the following design:

![](https://raw.githubusercontent.com/UC-Davis-molecular-computing/scadnano/master/doc-images/screenshot-initial.png)

The following Python script produces this design.

```python
import scadnano as sc
import modifications as mod


def main():
    # helices
    helices = [sc.Helix(max_offset=48), sc.Helix(max_offset=48)]

    # left staple
    stap_left_domain1 = sc.Domain(helix=1, forward=True, start=8, end=24)
    stap_left_domain0 = sc.Domain(helix=0, forward=False, start=8, end=24)
    stap_left = sc.Strand(domains=[stap_left_domain1, stap_left_domain0])

    # right staple
    stap_right_domain0 = sc.Domain(helix=0, forward=False, start=24, end=40)
    stap_right_domain1 = sc.Domain(helix=1, forward=True, start=24, end=40)
    stap_right = sc.Strand(domains=[stap_right_domain0, stap_right_domain1])
    stap_right.set_modification_5p(mod.biotin_5p)

    # scaffold
    scaf_domain1_left = sc.Domain(helix=1, forward=False, start=8, end=24)
    scaf_domain0 = sc.Domain(helix=0, forward=True, start=8, end=40)
    loopout = sc.Loopout(length=3)
    scaf_domain1_right = sc.Domain(helix=1, forward=False, start=24, end=40)
    scaf = sc.Strand(domains=[scaf_domain1_left, scaf_domain0, loopout, scaf_domain1_right], is_scaffold=True)

    # whole design
    design = sc.DNADesign(helices=helices, strands=[scaf, stap_left, stap_right], grid=sc.square)

    # deletions and insertions added to design are added to both strands on a helix
    design.add_deletion(helix=1, offset=20)
    design.add_insertion(helix=0, offset=14, length=1)
    design.add_insertion(helix=0, offset=26, length=2)

    # also assigns complement to strands other than scaf bound to it
    design.assign_dna(scaf, 'AACGT' * 18)

    return design


if __name__ == '__main__':
    design = main()
    design.write_scadnano_file(directory='output_designs')
```

Running the code above produces the `.dna` JSON file shown in the [web interface README](https://github.com/UC-Davis-molecular-computing/scadnano/blob/master/README.md#terminology). That section explains many of the terms used in the code (domain, helix, loopout, insertion, etc.).


## abbreviated syntax with chained methods
Instead of explicitly creating variables and objects representing each domain in each strand, there is a shorter syntax using chained method calls. Instead of the above, create only the helices first, then create the DNADesign. Then strands can be added using a shorter syntax, to describe how to draw the strand starting at the 5' end and moving to the 3' end. The following is a modified version of the above script using these chained methods

```python
import scadnano as sc
import modifications as mod


def main():
    # helices
    helices = [sc.Helix(max_offset=48), sc.Helix(max_offset=48)]

    # whole design
    design = sc.DNADesign(helices=helices, strands=[], grid=sc.square)

    # left staple
    design.strand(1, 8).to(24).cross(0).to(8)

    # right staple
    design.strand(0, 40).to(24).cross(1).to(40).with_modification_5p(mod.biotin_5p)

    # scaffold
    design.strand(1, 24).to(8).cross(0).to(40).loopout(1, 3).to(24).as_scaffold()

    # deletions and insertions added to design are added to both strands on a helix
    design.add_deletion(helix=1, offset=20)
    design.add_insertion(helix=0, offset=14, length=1)
    design.add_insertion(helix=0, offset=26, length=2)

    # also assigns complement to strands other than scaf bound to it
    design.assign_dna(design.scaffold, 'AACGT' * 18)

    return design


if __name__ == '__main__':
    design = main()
    design.write_scadnano_file(directory='output_designs')
```

Documentation is available in the [API docs](https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.DNADesign.strand).







## Documentation

Online documentation of the package API is located here:
https://scadnano-python-package.readthedocs.io





## Tutorial

A [tutorial](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/blob/master/tutorial/tutorial.md) shows how to create a "standard" 24-helix DNA origami rectangle using the scripting library.





## Other examples

Several example scripts are located in the 
[examples/](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/tree/master/examples) subfolder. 
Their output is contained in the 
[examples/output_designs/](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/tree/master/examples/output_designs) subfolder.
