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


function build_tests()
{
   local num_tests="${#sol_p_list[@]}"
   local test_id= sol_p_lst= ir_order_lst= exe_sfx_lst=" " exe_lst=
   for test_id; do
      set_test_params $test_id || continue
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
   $dry_run make -j $num_proc_build $test_name "${make_extra[@]}" || exit 1
}


function split_power2()
{
   local number="$2" outvar="$1" n1= n2= n3=
   [[ "$geom" = "Geometry::SQUARE" ]] && {
      (( n2 = 2**(number/2) ))
      (( n1 = n2*(2**(number%2>0?1:0)) ))
      eval "$outvar=($n1 $n2)"
      return 0
   }
   [[ "$geom" = "Geometry::CUBE" ]] && {
      (( n3 = 2**(number/3) ))
      (( n2 = n3*(2**(number%3>1?1:0)) ))
      (( n1 = n3*(2**(number%3>0?1:0)) ))
      eval "$outvar=($n1 $n2 $n3)"
      return 0
   }
   echo "Invalid Geometry."
   exit 1
}


function make_mesh_file()
{
   [[ ${#mesh_nxyz[@]} -eq 2 ]] && {
      printf -v mesh_file "quad-%03dx%03d.mesh" "${mesh_nxyz[@]}"
      [[ -n "$dry_run" ]] || [[ -e "$mesh_file" ]] || {
         printf "%s" "\
MFEM INLINE mesh v1.0

type = quad
nx = ${mesh_nxyz[0]}
ny = ${mesh_nxyz[1]}
sx = 1.0
sy = 1.0
"
      } > "$test_exe_dir/$mesh_file"
      return 0
   }
   [[ ${#mesh_nxyz[@]} -eq 3 ]] && {
      printf -v mesh_file "hex-%02dx%02dx%02d.mesh" "${mesh_nxyz[@]}"
      [[ -n "$dry_run" ]] || [[ -e "$mesh_file" ]] || {
         printf "%s" "\
MFEM INLINE mesh v1.0

type = hex
nx = ${mesh_nxyz[0]}
ny = ${mesh_nxyz[1]}
nz = ${mesh_nxyz[2]}
sx = 1.0
sy = 1.0
sz = 1.0
"
      } > "$test_exe_dir/$mesh_file"
      return 0
   }
   echo "Invalid mesh_nxyz variable."
   exit 1
}


function run_test()
{
   set_mpi_options
   local test_name_sfx="${test_name}${suffix}"
   local common_args="-no-vis $mesh_opt -rs $ser_ref -rp $par_ref -pc none"
   local num_args="${#args_list[@]}" args= all_args=()
   for ((i = 0; i < num_args; i++)) do
      args="${args_list[$i]}"
      all_args=($common_args $args "${test_extra_args[@]}")
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


function configure_tests()
{

# Set variables used by the functions 'build_tests', 'run_test',
# 'set_test_params', and 'run_tests_if_enabled':
test_name=mass
# problem: 0 - diffusion, 1 - mass
problem=${problem:-1}
dim=${dim:-2}
case "$dim" in
   2) geom="Geometry::SQUARE"
      mesh_s_list=( 17  17  18  17  16  16  15  14  17  17)
      ;;
   3) geom="Geometry::CUBE"
      mesh_s_list=( 15  15  16  15  14  13  13  12  15  15)
      ;;
   *) echo "Invalid dim = $dim. Stop."
      return 1
      ;;
esac
mesh_p=1
# test id:     0   1   2   3   4   5   6   7   8   9
sol_p_list=(   1   2   3   4   5   6   7   8   1   2)
ir_order_list=(0   0   0   0   0   0   0   0   3   5)
par_ref_list=( 2   1   0   0   0   0   0   0   2   1) # set below
enabled_tests="0   1   2   3   4   5   6   7   8   9"
# enabled_tests="2"
ser_ref=0
mesh_s_reduction_base=0     # used in run_tests_if_enabled()
mesh_s_reduction_limit=0    # used in run_tests_if_enabled()
mesh_s_max=29 # s=30 causes overflow in number of faces/edges
use_mpi_wtime=yes # leave empty for 'no'
rebuild_tests=no

local n=$num_proc_run
proc_t=0
while (( n >=2 )); do
   ((proc_t=proc_t+1))
   ((n=n/2))
done
(( num_proc_run == 2**proc_t )) || {
   echo "Invalid number of processors: $num_proc_run"; return 1; }
split_power2 proc_nxyz $proc_t
test_extra_args=(-c "${proc_nxyz[*]}")

case $num_nodes in
   1)
      par_ref_list=( 2   1   0   0   0   0   0   0   2   1)
      ;;
   8)
      par_ref_list=( 3   2   1   1   1   1   1   1   3   2)
      ;;
   64)
      par_ref_list=( 4   3   2   2   2   2   2   2   4   3)
      ;;
   512)
      par_ref_list=( 5   4   3   3   3   3   3   3   5   4)
      ;;
   *)
      printf "\nNumber of nodes = $num_nodes is not configured!\n"
      return 1
      ;;
esac

return 0

}

function set_test_params()
{
   [[ -n "$1" ]] && {
      sol_p="${sol_p_list[$1]}"
      ir_order="${ir_order_list[$1]}"
      mesh_s="${mesh_s_list[$1]}"
      (( mesh_s = mesh_s + mesh_s_shift ))
      par_ref="${par_ref_list[$1]}"
   }
   # echo "proc_t = $proc_t, mesh_s = $mesh_s, ser_ref = $ser_ref"
   if (( mesh_s + dim*(ser_ref+par_ref) < 0 )); then
      # invalid number of elements
      return 1
   fi
   if (( mesh_s + dim*(ser_ref+par_ref) > mesh_s_max )); then
      # too many elements
      return 1
   fi
   # adjust mesh_s as needed:
   if (( proc_t > mesh_s + dim*ser_ref )); then
      if (( proc_t <= mesh_s + dim*(ser_ref+par_ref) )); then
         while (( proc_t > mesh_s + dim*ser_ref )); do
            (( mesh_s = mesh_s + dim ))
            (( par_ref = par_ref - 1 ))
         done
      else
         # number of processors is larger than number of elements
         return 1
      fi
   fi
   suffix=_${geom#Geometry::}_p${sol_p}_m${mesh_p}
   (( ir_order != 0 )) && suffix+="_i${ir_order}"
   (( problem != 1 )) && suffix="_diff$suffix"
   split_power2 mesh_nxyz $mesh_s
   make_mesh_file
   mesh_opt="-m $mesh_file"
}

function run_tests_if_enabled()
{
   local test_id= enabled_tests_pat=" $enabled_tests " sft=
   for test_id; do
      if [[ -z "${enabled_tests_pat##* $test_id *}" ]]; then
         for ((sft = mesh_s_reduction_base;
               sft <= mesh_s_reduction_limit; sft++ )); do
            (( mesh_s_shift = -sft ))
            set_test_params $test_id && run_test || {
               if (( mesh_s+dim*(ser_ref+par_ref) >= 0 )); then
                  printf "Skipped test #$test_id with number of mesh elements"
                  echo " = $(( 2**(mesh_s+dim*(ser_ref+par_ref)) ))."
               fi
            }
         done
      fi
   done
}


function build_and_run_tests()
{

$dry_run cd "$test_dir"
configure_tests || return 1

build_tests $enabled_tests
echo

$dry_run cd "$test_exe_dir"
args_list=('-perf -mf')
total_memory_required_list=(8)  # guess-timates
run_tests_if_enabled 0 1 2 3 4 5 6 7 8 9

$dry_run make -f "$test_dir/makefile" clean-exec

}


test_required_packages="metis hypre mfem-patched"
mfem_patch_file="$test_dir/mass.patch"
