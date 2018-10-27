# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
# reserved. See files LICENSE and NOTICE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research is supported by the Exascale Computing Project
# (17-SC-20-SC), a collaborative effort of two U.S. Department of Energy
# organizations (Office of Science and the National Nuclear Security
# Administration) responsible for the planning and preparation of a capable
# exascale ecosystem, including software, applications, hardware, advanced
# system engineering and early testbed platforms, in support of the nation's
# exascale computing imperative.


function build_tests()
{
   local make_args="MFEM_DIR=${MFEM_DIR} BLD=${test_exe_dir}"
   quoted_echo make -f makefile-engines -j $num_proc_build $make_args
   [[ -n "$dry_run" ]] || make -f makefile-engines -j $num_proc_build $make_args
}


function run_test()
{
   set_mpi_options
   local common_args="-no-vis $mesh_opt -rs $ser_ref -rp $par_ref -o $sol_p"
   common_args+=" -p $problem"
   if [[ "$force_cuda_aware_mpi" == "1" ]]; then
      common_args+=" -cm"
   fi
   local all_args=($common_args "${test_extra_args[@]}")
   # ---=== TODO ===---
   total_memory_required="8"
   if [ -z "$dry_run" ]; then
      echo "Running test:"
      quoted_echo $mpi_run ./$test_name "${all_args[@]}"
      check_memory_req && {
         $mpi_run ./$test_name "${all_args[@]}" || \
         printf "\nError in the test, error code: $?\n\n"
      }
   else
      $dry_run $mpi_run ./$test_name "${all_args[@]}"
      check_memory_req
   fi
}


# Called one time by the main function, 'build_and_run_tests'
function configure_tests()
{

# Set variables used by the functions in this file.

force_cuda_aware_mpi=${force_cuda_aware_mpi:-0}

test_name=bp_engines
# problem: 0 - mass, 1 - diffusion
problem=${problem:-1}
# test id:     1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
sol_p_list=(   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)
ir_order_list=(0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0)
# enabled_tests_def="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16"
enabled_tests_def="1 2 3 4 5 6 7 8"
# 'enabled_tests' can be set on the 'go.sh' command line
enabled_tests="${enabled_tests:-$enabled_tests_def}"
mesh_file="${MFEM_DIR}/data/star.mesh"
mesh_dim=2
mesh_num_elem=20
ser_ref=0
max_dofs=${max_dofs:-3000000}
rebuild_tests=no

return 0

}


# Called by 'run_tests'
function set_test_params()
{
   local tid=$(($1-1))
   sol_p=${sol_p_list[$tid]}
   ir_order=${ir_order_list[$tid]}
   suffix=_${geom#Geometry::}_p${sol_p}_m${mesh_p}
   (( ir_order != 0 )) && suffix+="_i${ir_order}"
   mesh_opt="-m $mesh_file"
}


# Called one time by the main function, 'build_and_run_tests'
function run_tests()
{
   for test_id; do
      set_test_params $test_id || continue
      local dofs_0=$((mesh_num_elem*sol_p**mesh_dim))
      local ref_fac=$((2**mesh_dim))
      for ((par_ref = 0; dofs_0*ref_fac**par_ref < max_dofs; par_ref++)); do
         run_test
      done
   done
}


# Main function
function build_and_run_tests()
{

$dry_run cd "$test_dir"
configure_tests || return 1

echo "Enabled tests: $enabled_tests"
build_tests $enabled_tests || return 1
echo

[[ -n "$build_only" ]] && return

$dry_run cd "$test_exe_dir"
run_tests $enabled_tests

}


# test_required_packages="metis hypre occa-fork mfem-engines"
test_required_packages="metis hypre occa mfem-engines"
