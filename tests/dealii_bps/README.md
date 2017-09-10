
This directory contains scripts for setting up, running and postprocessing
the [deal.II CEED BPs](https://github.com/kronbichler/ceed_benchmarks_dealii)
written by Martin Kronbichler.

## Running the scripts

Example:
```sh
../../go.sh -c linux -m gcc -r bp1.sh -n 16 -p 16
```
where `-n 16` is the total number of processors and `-p 16` is the number of
processors per node.

Multiple processor configuration can be run with:
```sh
../../go.sh -c linux -m gcc -r bp1.sh -n "16 32 64" -p "16 32 64"
```

## Post-processing the results

First, save the output of the run to a file:
```sh
../../go.sh -c linux -m gcc -r bp1.sh -n 16 -p 16 > run_001.txt
```
and then use one of the `postprocess-plot-*.py` scripts (which require
the python package matplotlib), e.g.:
```sh
python postprocess-plot-2.py run_001.txt
```
Note that the `postprocess-*.py` scripts can read multiple files at a time just
by listing them on the command line and also read the standard input if no files
were specified on the command line.

### Plotting performance comparisons

The script `postprocess-plot-4.py` can be used to plot the performance ratios
between two sets of data, e.g. vector (BP2) vs scalar (BP1), or two compilers:
XLC vs GCC. Sample usage:
```sh
python postprocess-plot-4.py bp1_gcc_n16.txt bp2_gcc_n16.txt
python postprocess-plot-4.py bp1_gcc_n16.txt bp1_xlc_n16.txt
```
