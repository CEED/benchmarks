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

test_name=bp_main
# problem: 0 - diffusion, 1 - mass
problem=${problem:-0}
geom=Geometry::SQUARE
mesh_p=1
# quadrature type: 0 - Gauss, 1 - Gauss-Lobatto
ir_type=${ir_type:-0}
vdim=${vdim:-1}
vec_layout=${vec_layout:-}
# test id:     0   1   2   3   4   5   6   7
sol_p_list=(   1   2   3   4   5   6   7   8)
ir_order_list=(3   5   7   9  11  13  15  17)
enabled_tests_def="0 1 2 3 4 5 6 7"
(( ir_type != 0 )) && enabled_tests_def="0 1 2 3 4 5 6 7"
# enabled_tests_def="0"
enabled_tests="${enabled_tests:-$enabled_tests_def}"
mesh_file="${MFEM_DIR}/data/star.mesh"
mesh_dim=2
mesh_num_elem=20
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
   suffix=_${geom#Geometry::}_p${sol_p}_m${mesh_p}
   (( ir_order != 0 )) && suffix+="_i${ir_order}"
   mesh_opt="-m $mesh_file"
}


function build_tests()
{
   local test_id= sol_p_lst= ir_order_lst= exe_sfx_lst=" " exe_lst=
   for test_id; do
      set_test_params $test_id

      # check if suffix is already in the list
      if [[ -n "${exe_sfx_lst##* $suffix *}" ]]; then
         sol_p_lst+=" $sol_p"
         ir_order_lst+=" $ir_order"
         exe_sfx_lst+="$suffix "
         exe_lst+=" ${test_name}$suffix"
      fi
   done
   exe_sfx_lst="${exe_sfx_lst% }"
   [[ -n "$sol_p_lst" ]] || return 0
   case "$rebuild_tests" in
      yes|Yes|YES)
         $dry_run cd "$test_exe_dir"
         $dry_run make -f "$test_dir/makefile" clean
         $dry_run rm -f ${exe_lst:1}
         $dry_run cd "$test_dir"
         ;;
   esac
   local make_extra=("geom=$geom" "mesh_p=$mesh_p" "sol_p=${sol_p_lst:1}")
   make_extra=("${make_extra[@]}" "problem=$problem")
   make_extra=("${make_extra[@]}" "use_mpi_wtime=$use_mpi_wtime")
   make_extra=("${make_extra[@]}" "ir_order=${ir_order_lst:1}")
   make_extra=("${make_extra[@]}" "exe_suffix=${exe_sfx_lst:1}")
   make_extra=("${make_extra[@]}" "EXTRA_CXXFLAGS=$TEST_EXTRA_CFLAGS")
   make_extra=("${make_extra[@]}" "MFEM_DIR=$MFEM_DIR")
   make_extra=("${make_extra[@]}" "BLD=$test_exe_dir/")
   make_extra=("${make_extra[@]}" "EXTRA_INCFLAGS=-I$MFEM_SOURCE_DIR")
   quoted_echo make -j $num_proc_build $test_name "${make_extra[@]}"
   [[ -n "$dry_run" ]] || make -j $num_proc_build $test_name "${make_extra[@]}"
}


function run_test()
{
   set_mpi_options
   local test_name_sfx="${test_name}${suffix}"
   local common_args="-no-vis $mesh_opt -rs $ser_ref -rp $par_ref -pc none"
   local num_args="${#args_list[@]}" args= all_args=()
   for ((i = 0; i < num_args; i++)) do
      args="${args_list[$i]}"
      all_args=($common_args $args)
      total_memory_required="${total_memory_required_list[$i]}"
      if [ -z "$dry_run" ]; then
         echo "Running test:"
         quoted_echo $mpi_run ./$test_name_sfx "${all_args[@]}"
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

echo "Enabled tests: $enabled_tests"
build_tests $enabled_tests || return 1
echo

[[ -n "$build_only" ]] && return

$dry_run cd "$test_exe_dir"
args_list=('-perf -mf')
total_memory_required_list=(8)  # guess-timates
run_tests $enabled_tests

$dry_run make -f "$test_dir/makefile" clean-exec

}


test_required_packages="metis hypre mfem"
