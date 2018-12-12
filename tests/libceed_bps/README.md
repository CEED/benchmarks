## Running the scripts

Example:
```sh
../../go.sh -c linux -m gcc -r bp1.sh -n 16 -p 16
```
where `-n 16` is the total number of nodes and `-p 16` is the number of
processors per node.

Multiple processor configuration can be run with:
```sh
../../go.sh -c linux -m gcc -r bp1.sh -n "16 32 64" -p "16 32 64"
```

The following variables can be set on the command line:
* `ceed=<libceed-device-spec>`, e.g. `ceed=/cpu/self/ref`; the default value is
  `/cpu/self/blocked`.
* `max_dofs_node=<number>`, e.g. `max_dofs_node=1000000` - this sets the upper
  bound of the problem sizes, per compute node; the default value is 3*2^20.
* `max_p=<number>`, e.g. `max_p=12` - this sets the highest degree for which the
  tests will be run (the lowest degree is 1); the default value is 8.
* `build_only=<string>`, e.g. `build_only=1` - if the string is not empty, the
  execution will stop after the executables are built; the default is the empty
  string. This option is useful for pre-building the required packages and the
  executables before running the actual tests.
* `LIBCEED_DIR=<directory>`, e.g. `LIBCEED_DIR=$HOME/code/libceed` - if this
  variable is set, the test will be linked against the given libCEED library,
  instead of the libCEED library that will be otherwise built by the script.

To remove a previous build of the test executable, e.g. to force a re-build,
simply delete the directory `../../builds/<config>_<compiler>/libceed_bps`.

## Post-processing the results

First, save the output of the run to a file:
```sh
../../go.sh -c linux -m gcc -r bp1.sh -n 16 -p 16 > run_001.txt
```
and then use the `postprocess-plot-2.py` script (which requires the python
package matplotlib) or the `postprocess-table.py` script, e.g.:
```sh
python postprocess-plot-2.py run_001.txt
```
The plot ranges and some other options can be adjusted by editting the values
in the beginning of the script `postprocess-plot-2.py`.

Note that the `postprocess-*.py` scripts can read multiple files at a time just
by listing them on the command line and also read the standard input if no files
were specified on the command line.
