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


function configure_tests()
{

test_name=laghos
# problem: 0 - diffusion, 1 - mass
problem=${problem:-1}
geom=Geometry::SQUARE
mesh_p=1
# quadrature type: 0 - Gauss, 1 - Gauss-Lobatto
ir_type=${ir_type:-0}
vdim=${vdim:-1}
vec_layout=${vec_layout:-}
# test id:     0   1   2
sol_p_list=(   2   3   4)
ir_order_list=(5   7   9)
enabled_tests_def="0 1 2"
(( ir_type != 0 )) && enabled_tests_def="0 1 2"
# enabled_tests_def="0"
enabled_tests="${enabled_tests:-$enabled_tests_def}"
mesh_file="../laghos/data/square01_quad.mesh"
mesh_dim=2
mesh_num_elem=2
ser_ref=0
max_dofs=${max_dofs:-3000000}
use_mpi_wtime=yes # leave empty for 'no'
rebuild_tests=no

return 0

}


function set_test_params()
{
   sol_p=${sol_p_list[$1]}
   ir_order=${ir_order_list[$1]}
   mesh_opt="-m $mesh_file"
}


function run_test()
{
   set_mpi_options
   local test_name_sfx="../laghos/${test_name}"
   local common_args="-p 1 -ok $sol_p -ot $((sol_p-1)) -s 7 -no-vis -cfl 1e-07 -cgt 0 -cgm 15 -ms 1 $mesh_opt -rs $ser_ref -rp $par_ref"
   local num_args="${#args_list[@]}" args= all_args=()
   for ((i = 0; i < num_args; i++)) do
      args="${args_list[$i]}"
      all_args=($common_args $args)
      total_memory_required="${total_memory_required_list[$i]}"
      if [ -z "$dry_run" ]; then
         echo "Running test:"
         quoted_echo $mpi_run $test_name_sfx "${all_args[@]}"
         check_memory_req && {
            $mpi_run ./$test_name_sfx "${all_args[@]}" || \
            printf "\nError in the test, error code: $?\n\n"
         }
      else
         $dry_run $mpi_run ./$test_name_sfx "${all_args[@]}"
         check_memory_req
      fi
   done
}


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


function build_and_run_tests()
{

$dry_run cd "$test_dir"
configure_tests || return 1

[[ -n "$build_only" ]] && return

$dry_run cd "$test_exe_dir"
args_list=('')
total_memory_required_list=(8)  # guess-timates
run_tests $enabled_tests

$dry_run make -f "$test_dir/makefile" clean-exec

}


test_required_packages="metis hypre mfem laghos"
