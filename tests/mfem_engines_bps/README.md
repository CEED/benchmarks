## Running the scripts

Example:
```sh
../../go.sh -c sierra -m xlc -r bp_main.sh -n 16 -p 4
```
where `-n 16` is the total number of MPI tasks and `-p 4` is the number of
tasks per node.

Multiple processor configuration can be run with:
```sh
../../go.sh -c sierra -m xlc -r bp_main.sh -n "16 32 64" -p "4 4 4"
```

## Post-processing the results

First, save the output of the run to a file:
```sh
../../go.sh -c sierra -m xlc -r bp_main.sh -n 16 -p 4 > run_001.txt
```
and then use one of the `postprocess-plot-*.py` scripts (which require
the python package matplotlib) or the `postprocess-table.py` script, e.g.:
```sh
python postprocess-plot-2.py run_001.txt
```
Note that the `postprocess-*.py` scripts can read multiple files at a time just
by listing them on the command line and also read the standard input if no files
were specified on the command line.
