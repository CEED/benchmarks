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

TODO

### Plotting performance comparisons

TODO
