## Running the scripts

Example:
```sh
../../go.sh -c linux -m gcc -r bp1_v1.sh -n 16 -p 16
```
where `-n 16` is the total number of processors and `-p 16` is the number of
processors per node. Note that `bp1_v1.sh` (and `bp3_v1.sh`, which uses
`bp1_v1.sh`) requires a power of 2 for the number of processors.

Multiple processor configuration can be run with:
```sh
../../go.sh -c linux -m gcc -r bp1_v1.sh -n "16 32 64" -p "16 32 64"
```

## Post-processing the results

First, save the output of the run to a file:
```sh
../../go.sh -c linux -m gcc -r bp1_v1.sh -n 16 -p 16 > run_001.txt
```
and then use one of the `postprocess-plot-*.py` scripts (which require
the python package matplotlib) or the `postprocess-table.py` script, e.g.:
```sh
python postprocess-plot-1.py run_001.txt
```
Note that the `postprocess-*.py` scripts can read multiple files at a time just
by listing them on the command line and also read the standard input if no files
were specified on the command line.
