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

function copy_mesh()
{
    $dry_run cp $test_dir/inline-hex.mesh $test_exe_dir/
}

function setup_occa_test()
{
    # Set test_name before calling this function
    [ -z $test_name ] && { echo "test_name is unset" && return; }

    # Export OCCA_CACHE_DIR
    export OCCA_CACHE_DIR=${OCCA_CACHE_DIR:-$test_exe_dir/.occa}

    # Build the executable
    $dry_run cd "$test_dir"

    $dry_run make "$test_exe_dir/$test_name" "MFEM_DIR=$MFEM_DIR" "BLD=$test_exe_dir/"
    [ -z $dry_run ] && {
        occa info && occa clear -ay && make cache-kernels "MFEM_DIR=$MFEM_DIR" "MFEM_SOURCE_DIR=$MFEM_SOURCE_DIR"
    }

    # Copy the mesh over
    copy_mesh
}

function setup_baseline_test()
{
    # Set test_name before calling this function
    [ -z $test_name ] && { echo "test_name is unset" && return; }

    # Build the executable
    $dry_run cd "$test_dir"

    $dry_run make "$test_exe_dir/$test_name" "MFEM_DIR=$MFEM_DIR" "BLD=$test_exe_dir/"

    # Copy the mesh over
    copy_mesh
}
