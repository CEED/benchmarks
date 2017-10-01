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

source tests/mfem_pao_bps/bps.sh

function build_and_run_tests()
{
    # Copy the mesh into place
    [ ! -r $test_exe_dir/inline-hex.mesh ] &&
        $dry_run cp $test_dir/inline-hex.mesh $test_exe_dir/
    mesh_dim=3

    # Build the executable
    $dry_run cd $test_exe_dir &&
        $dry_run make $test_exe_dir/bps MFEM_DIR=$MFEM_DIR BLD=$test_exe_dir/

    # Constant parameters
    quadorder=-1
    basistype=G
    maxiter=100
    print_level=1
    write_solution=0

    # Fixed Parameters
    bp=1
    mf=0

    # Loop through cases
    if (( bp % 2 )); then vdim=1; else vdim=$mesh_dim; fi
    for num_proc_run in 1 2 5 10 20; do
        for order in 1 3 5 7; do
            for ser_ref in {0..10}; do
                for par_ref in 0; do
                    (( dofs = 2**(mesh_dim + 3*(ser_ref+par_ref)) * (order**mesh_dim) * vdim ))
                    if (( dofs > 5000 )) && (( dofs < 1000000 )); then
                        run_test
                    fi
                done
            done
        done
    done
}

test_required_packages="metis hypre mfem"
mfem_git_branch="pa-oper-dev"
