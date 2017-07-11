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

# function configure_tests()
# {
#    problem=1
#    sol_p_list=(1 2 3 4 5 6 7 8 9)
   
# }

# function build_tests()
# {
# }

function build_and_run_tests()
{

   # Parameters
   local problem=4
   local sol_p=3
   local pc=lor #lumpedmass

   # Build test
   $dry_run cd "$test_dir"
   local make_extra=("PROBLEM=$problem" "SOL_P=$sol_p" "USE_MPI_WTIME=1")
   make_extra=("${make_extra[@]}" "EXTRA_CXXFLAGS=$TEST_EXTRA_CFLAGS")
   make_extra=("${make_extra[@]}" "MFEM_DIR=$MFEM_DIR")
   make_extra=("${make_extra[@]}" "BLD=$test_exe_dir/")
   quoted_echo make -j $num_proc_build "${make_extra[@]}"
   [[ -n "$dry_run" ]] || make -j $num_proc_build "${make_extra[@]}"

   # Run test
   $dry_run cd "$test_exe_dir"
   set_mpi_options
   local test_name=bp${problem}_solp${sol_p}
   local all_args="--no-visualization"
   all_args="$all_args -e 36"
   all_args="$all_args --preconditioner $pc"
   if [ -z "$dry_run" ]; then
      echo "Running test:"
      quoted_echo $mpi_run ./$test_name $all_args
      $mpi_run ./$test_name $all_args
   fi

   # Clean up test
   $dry_run make -f "$test_dir/makefile" clean-exec "${make_extra[@]}"

}


test_required_packages="metis hypre mfem"
