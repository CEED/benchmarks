## Bake-off/benchmark problems and other tests with Nek5000

### Required tools

* bash
* make
* git
* wget or curl

### Sample use

First, you may need to create a configuration for your machine, e.g. by copying
and modifying one of the existing `machine-configs/*.sh` files.

```sh
./go.sh --config vulcan --compiler gcc --build "nek5000"
./go.sh --config vulcan --compiler gcc --run tests/nek5000_bps/bp1/bp1.sh
```

To see a list of the available configs use `./go.sh` or generally use
`./go.sh --help` for help. These configs correspond to the scripts
`machine-configs/<name>.sh`.

To see the available compilers for a config use `./go.sh --config <name>`.
