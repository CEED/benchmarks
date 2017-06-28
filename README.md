# CEED Benchmarks

This repository contains bake-off/benchmark problems and other tests for
performance evaluation of high-order kernels on HPC architectures developed in
the ECP co-design [Center for Efficient Exascale Discretizations
(CEED)](http://ceed.exascaleproject.org).

The repository is part of the [CEED software
suite](http://ceed.exascaleproject.org/software/), a collection of software
benchmarks, miniapps, libraries and APIs for efficient exascale discretizations
based on high-order finite element and spectral element methods.  See
http://github.com/ceed for more information and source code availability.

The CEED research is supported by the [Exascale Computing
Project](https://exascaleproject.org/exascale-computing-project) (17-SC-20-SC),
a collaborative effort of two U.S. Department of Energy organizations (Office of
Science and the National Nuclear Security Administration) responsible for the
planning and preparation of a [capable exascale
ecosystem](https://exascaleproject.org/what-is-exascale/), including software,
applications, hardware, advanced system engineering and early testbed platforms,
in support of the nation’s [exascale computing
imperative](https://obamawhitehouse.archives.gov/the-press-office/2015/07/29/executive-order-creating-national-strategic-computing-initiative).

For more details on the CEED benchmarks see http://ceed.exascaleproject.org/bps/

## Building and Running the Benchmarks

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

Equivalent short version:

```sh
./go.sh -c vulcan -m gcc -b "metis hypre mfem"
./go.sh -c vulcan -m gcc -r tests/mfem_bps/bp1_v1.sh
```

To see a list of the available configs use `./go.sh` or generally use
`./go.sh --help` for help. These configs correspond to the scripts
`machine-configs/<name>.sh`.

To see the available compilers for a config use `./go.sh --config <name>`.

## Copyright

The following copyright applies to each file in the CEED software suite, unless
otherwise stated in the file:

Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at the
Lawrence Livermore National Laboratory. LLNL-CODE-XXXXXX. All Rights reserved.

