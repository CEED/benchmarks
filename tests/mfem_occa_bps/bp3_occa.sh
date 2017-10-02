# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-XXXXXX. All Rights
# reserved. See file LICENSE for details.
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

if [[ -z "$root_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi

source "$test_dir/setup.sh.inl"

function build_and_run_tests()
{
    # This test name
    test_name=bp3

    # Setup and build the executable
    setup_occa_test

    # The test definition
    orders=(1 2 3 4 5 6 7 8)
    ref_min=(1 1 1 1 0 0 0 0)
    ref_max=(6 5 4 4 4 3 3 3)

    mode="mode:'Serial'"
    case $OCCA_MODE in
        Serial) mode="mode:'Serial'";;
        CUDA)   mode="mode:'CUDA',deviceID:0";;
        OpenMP) mode="mode:'OpenMP'";;
    esac

    # Run the tests
    $dry_run cd "$test_exe_dir"
    for ((i=0; i<"${#orders[@]}"; i++)); do
        args="-d $mode --order ${orders[$i]} --no-acro --mesh inline-hex.mesh"
        $dry_run $CLI_PREFIX ./$test_name $args -i 5 -r 1 -ov &>/dev/null
        if (( ${orders[$i]} < 3 )); then
            $dry_run $CLI_PREFIX ./$test_name $args -i 5 -r 1 -ov -q -2 &>/dev/null
        fi
        for ((ref=${ref_min[$i]}; ref<=${ref_max[$i]}; ref++)); do
            $dry_run $CLI_PREFIX ./$test_name $args --ref-levels $ref --occa-verbose
            if (( ${orders[$i]} < 3 )); then
                $dry_run $CLI_PREFIX ./$test_name $args --ref-levels $ref --occa-verbose --quad-add -2
            fi
        done
    done
}


test_required_packages="metis hypre occa mfem-occa"
