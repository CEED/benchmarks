# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
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

function configure_tests()
{
   problem=${problem:-1}
   sol_p_list=(   1  2  3   4   5   6   7   8  1  2)
   ir_order_list=(5  7  9  11  13  15  17  19  3  5)
   el_per_proc_list=(1 2 4 8 12 18 27 36)
   pc=${pc:-none}
}

function build_tests()
{

   # Move to test directory
   $dry_run cd "$test_dir"

   # Initialize build arguments
   local make_extra=("PROBLEM=${problem}" "sol_p=${sol_p_list[*]}" "ir_order=${ir_order_list[*]}" "USE_MPI_WTIME=1")
   make_extra=("${make_extra[@]}" "EXTRA_CXXFLAGS=$TEST_EXTRA_CFLAGS")
   make_extra=("${make_extra[@]}" "MFEM_DIR=$MFEM_DIR")
   make_extra=("${make_extra[@]}" "BLD=$test_exe_dir/")

   # Clear previous builds if needed
   case "$rebuild_tests" in
      yes|Yes|YES)
         $dry_run make "${make_extra[@]}" clean
   esac

   # Build tests
   quoted_echo make -j $num_proc_build "${make_extra[@]}" all
   [[ -n "$dry_run" ]] || make -j $num_proc_build "${make_extra[@]}" all

}

function run_tests()
{

   # Move to test executable directory
   $dry_run cd "$test_exe_dir"

   # Initialize MPI options
   set_mpi_options

   # Iterate through test configurations
   local num_exes=${#sol_p_list[@]}
   local el_per_proc_list_size=${#el_per_proc_list[@]}
   for ((j = 0; j < el_per_proc_list_size; j++)) do
      for ((i = 0; i < num_exes; i++)) do
         local test_name=bp${problem}_solp${sol_p_list[i]}_irorder${ir_order_list[i]}
         local all_args="--no-visualization"
         all_args="$all_args --num-el-per-proc ${el_per_proc_list[j]}"
         all_args="$all_args --preconditioner $pc"
         if [ -z "$dry_run" ]; then
            echo "Running test:"
            quoted_echo $mpi_run ./$test_name $all_args
            $mpi_run ./$test_name $all_args
         fi
      done
   done
    
}

function build_and_run_tests()
{
   configure_tests || return 1
   build_tests || return 1
   [[ -n "$build_only" ]] && return
   run_tests || return 1
}

test_required_packages="metis hypre mfem"
