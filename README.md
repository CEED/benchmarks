## Bake-off/benchmark problems and other tests

### Required tools

* bash
* make, cmake (for some packages)
* git
* wget or curl

### Sample use

First, you may need to create a configuration for your machine, e.g. by copying
and modifying one of the existing `machine-configs/*.sh` files.

```sh
./go.sh --config vulcan --compiler gcc --build "metis hypre mfem"
./go.sh --config vulcan --compiler gcc --run tests/mfem_bps/bp1_v1.sh
```

To see a list of the available configs use `./go.sh` or generally use
`./go.sh --help` for help. These configs correspond to the scripts
`machine-configs/<name>.sh`.

To see the available compilers for a config use `./go.sh --config <name>`.
