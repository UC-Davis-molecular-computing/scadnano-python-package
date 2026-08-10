# AutoStaple

These are all the instructions I followed to get auto-staple working for GRID-type DNA. It is recommended you read the `README.md` for [cadnano2](https://github.com/douglaslab/cadnano2/blob/main/README.md) as well as follow the tutorial for the [scadnano-python-package](/tutorial/tutorial.md).

## Setup

This setup will be split into two parts -- one for [scadnano-python-package](https://github.com/UC-Davis-molecular-computing/scadnano-python-package) and [cadnano2](https://github.com/douglaslab/cadnano2).

This assumes you have Python 3.14, or higher, installed.

You can use a Python debugger of your choice, but I recommend to use `Pycharm`. It is a powerful IDE that is free with a `.edu` email address through [GitHub Education](https://github.com/education).

### Setup - cadnano2

In order to reverse-engineer the *AutoStaple* algorithm from `cadnano2` into `scadnano-python-package` we need to first setup `cadnano2` on our local machines.

1. Clone the [cadnano2](https://github.com/douglaslab/cadnano2) repository.

```shell
git clone https://github.com/douglaslab/cadnano2
cd cadnano2
```

2. Follow the setup instructions in the [README.md](https://github.com/douglaslab/cadnano2?tab=readme-ov-file#development) file in the `cadnano2`'s repository root for setting up the development environment for your machine.

3. Create a file in the root of the `cadnano2` repository titled `debug_launcher.py` with the following contents.

```python
from cadnano2.main import main

if __name__ == '__main__':
    main()

# Uncomment this if you will run the debug_launcher.py file directly.
# main()
```

This will allow us to run `cadnano2` in the debugger. Otherwise, you would need to recompile cadnano2 each time you made a change, and be unable to run through it line-by-line.

### Setup - scadnano-python-package

> Ensure that this is done in a separate terminal tab/window from the cadnano2 steps above.

1. Clone the repository.

```shell
git clone https://github.com/UC-Davis-molecular-computing/scadnano-python-package
cd scadnano-python-package
```

2. It is recommended to create a [virtual environment](https://docs.python.org/3/library/venv.html) of your choice for python. It is the same way as creating a virtual environment for `cadnano2`. You can also follow the instructions below.

```shell
pip3 install venv
python3 -m venv env
```

* MacOS/Linux

```shell
source env/bin/activate
```

* Windows

```shell
./env/Scripts/Activate.ps1
```


3. Run the following command to install dependencies.

```shell
pip install -e .
```

4. Inside the repository root, create a file called `test_autostaple.py` file within the `scadnano` directory. Use the contents below as a starting point.

```python
import scadnano as sc

if __name__ == '__main__':
    # Load an existing .sc design file, relative to the test_autostaple.py file.
    design1 = sc.Design.from_scadnano_file(filename='input/helix_test.sc')
    design2 = sc.Design.from_scadnano_file(filename='input/helix_test_2.sc')
    design3 = sc.Design.from_scadnano_file(filename='input/honeycomb_test.sc')

    # Run autostaple on it
    design1.autostaple()
    design2.autostaple()
    design3.autostaple()

    # Save the result, relative to teh test_autostaple.py file
    design1.write_scadnano_file(filename='helix_result.sc', directory='output')
    design2.write_scadnano_file(filename='helix_result_2.sc', directory='output')
    design3.write_scadnano_file(filename='honeycomb_result.sc', directory='output')
```

This will allow us to use the debugger with scadnano.

## Reverse Engineering AutoStaple

What helped me the most was wrapping my head around the differences between cadnano2 and scadnano. It seems like what is commonly known as a `Strand` in cadnano2 is a `Domain` in scadnano. (Domains in scadnano HAVE strands.)

Furthermore, it is important to note that `cadnano2` does things through "commands" whereas the implementation I have in `scadnano-python-package` is done directly on the part.

Initially, I used default designs from [scadnano.org](https://scadnano.org) using `File -> Load Example` as a basis.

For example, if I wanted to test the AutoStaple implementation for a basic square grid, I would do the following.

### Workflow

1. Head to [scadnano.org](https://scadnano.org).

2. Load an Example File: `File -> Load Example -> 6_helix_origami_rectangle`.

3. Export a scadnano file: `File -> Save`.

4. Export a cadnano2 file: `Export -> cadnano v2`.

5. Run cadnano2 in the `debugger` with breakpoints placed mirroring the breakpoints in scadnano-python-package.

6. Step through the code, confirming that the variables and algorithms mirror cadnano2.

Code for AutoStaple is located in the `autostaple` method in `scadnano.py` and `part.py` for scadnano-python-package and cadnano2 respectively.

## Work So Far

So far, I have transferred the majority of the AutoStaple code from cadnano2 to scadnano-python-package.

AutoStaple for square patterns seems to be fully working. However, AutoStaple for honeycomb patterns seems to be partially working.

Most of the "helper" functions I added start with an underscore (`_`). They have the same names as the functions in cadnano2.

### Areas of Interest

If you take a look at the `_getVirtualHelixNeighbors` in scadnano-python-package, and compare it to `getVirtualHelixNeighbors` in cadnano2, I have not yet implemented the separate rows. I suspect this is what is preventing AutoStaple from working properly for honeycomb grids.
