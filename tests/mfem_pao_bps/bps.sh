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

function run_test()
{
    set_mpi_options

    test_name=bps

    all_args=(--bakeoff $bp \
                     --refine-serial $ser_ref \
                     --refine-parallel $par_ref \
                     --order $order \
                     --basis-type $basistype \
                     --quadrature-order $quadorder \
                     --max-iter $maxiter \
                     --print-level $print_level)
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

function build_and_run_tests()
{
    # Copy the mesh into place
    [ ! -r $test_exe_dir/inline-hex.mesh ] && $dry_run cp $test_dir/inline-hex.mesh $test_exe_dir/
    mesh_dim=3

    # Build the executable
    $dry_run cd $test_exe_dir && $dry_run make $test_exe_dir/bps MFEM_DIR=$MFEM_DIR BLD=$test_exe_dir/


    # Constant parameters
    basistype=G
    maxiter=100
    print_level=1
    write_solution=0
    quadorder=-1

    # Loop through cases
    bp=1
    mf=0

    if (( bp % 2 )); then vdim=1; else vdim=$mesh_dim; fi
    for num_proc_run in 1; do
        for order in 1 3; do
            for ser_ref in {0..10}; do
                for par_ref in 0; do
                    (( dofs = 2**(mesh_dim + 3*(ser_ref+par_ref)) * (order**mesh_dim) * vdim ))
                    if (( dofs > 5 )) && (( dofs < 100000 )); then
                        run_test
                    fi
                done
            done
        done
    done
}

test_required_packages="metis hypre mfem"
mfem_git_branch="pa-oper-dev"
