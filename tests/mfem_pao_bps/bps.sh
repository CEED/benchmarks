# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
# reserved. See files LICENSE and NOTICE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research is supported by the Exascale Computing Project (17-SC-20-SC)
# a collaborative effort of two U.S. Department of Energy organizations (Office
# of Science and the National Nuclear Security Administration) responsible for
# the planning and preparation of a capable exascale ecosystem, including
# software, applications, hardware, advanced system engineering and early
# testbed platforms, in support of the nation's exascale computing imperative.

# This script by itself does not define a test case. It is included in
# the other bp*.sh scripts in this directory.

function run_test()
{
    local test_name=bps

    all_args=( --bakeoff $bp \
               --refine-serial $ser_ref \
               --refine-parallel $par_ref \
               --order $order \
               --basis-type $basistype \
               --quadrature-order $quadorder \
               --max-iter $maxiter \
               --print-level $print_level )

    extra_args=()
    if (( $mf )); then
        extra_args+=(--matrix-free)
    else
        extra_args+=(--assembly)
    fi

    if (( $write_solution )); then
        extra_args+=(--write-solution)
    else
        extra_args+=(--no-write-solution)
    fi

    quoted_echo $mpi_run ./$test_name "${all_args[@]}" "${extra_args[@]}"
    $dry_run $mpi_run ./$test_name "${all_args[@]}" "${extra_args[@]}"
}
