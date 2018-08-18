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

*TODO*
