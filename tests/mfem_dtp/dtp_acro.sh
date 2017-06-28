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

source "$test_dir/dtp.sh.inl"


function build_and_run_tests()
{
    # Run the tests
    local test_name=dtp_occa

    setup_occa_test

    $dry_run cd "$test_exe_dir"
    for i in "${!case_orders[@]}"; do
        $dry_run ./$test_name -d "mode:'CUDA',deviceID:0" --order "${case_orders[$i]}" --ref-levels "${case_refs[$i]}" -ac --preconditioner none --mesh "fichera.mesh"
    done
}


test_required_packages="metis hypre occa acrotensor mfem-occa"
