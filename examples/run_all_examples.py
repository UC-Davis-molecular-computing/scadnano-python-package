import os
import sys
import subprocess
import importlib
import scadnano as sc


def main() -> None:
    """
    Looks for .py files in the current directory. If they have a function called create_design(),
    it calls that to produce a design and write its output to directory output_designs
    using method Design.write_scadnano_file(); otherwise it runs the Python file as a
    command-line application using the subprocess module.
    """
    for filename in os.listdir("."):
        if filename.endswith(".py") and not sys.argv[0].endswith(filename):
            print(f"running {filename}")
            run(filename)


def run(filename: str) -> None:
    if not filename.endswith('.py'):
        print(f"  filename {filename} does not end with '.py'; skipping")
        return
    modulename = filename[:-3]
    print(f'  importing module {modulename}')
    module = importlib.import_module(modulename)
    if hasattr(module, "create_design"):
        print(f"  found create_design function in module {modulename}; running it and writing its output")
        design: sc.Design = module.create_design()
        design.write_scadnano_file(directory='output_designs',
                                   filename=modulename + f'.{sc.default_scadnano_file_extension}')
    else:
        print(f"  found no main function in module {modulename}; running as subprocess instead")
        try:
            retcode = subprocess.call(f"python {filename}", shell=True)
            if retcode < 0:
                print("Child was terminated by signal", -retcode, file=sys.stderr)
        except OSError as e:
            print("Execution failed:", e, file=sys.stderr)


if __name__ == "__main__":
    main()
