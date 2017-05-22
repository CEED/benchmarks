## MFEM+OCCA Examples

### Sample runs

```sh
../../go.sh -c mac -m clang -r ex1.sh
```

Edit the script `ex1.sh` to change the command line options for the example.

The same script can be used to build and run any of the MFEM+OCCA examples.

Note: the executable and any output will be placed in
`../../builds/mac_clang/mfem-occa/examples/occa`. Running the executable
directly may fail due to the different environment set by the chosen config,
compiler, package, and test scripts.
