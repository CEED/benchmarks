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


# Call graph of the functions in this file:
#
# +--+  build_and_run_tests
#    |
#    +--+  configure_tests
#    |  |
#    |  +--+  split3_power2
#    |
#    +--+  build_tests
#    |
#    +--+  run_tests
#       |
#       +--+  set_test_params
#       |  |
#       |  +--+  split3_power2
#       |  |
#       |  +--+  make_mesh_file
#       |
#       +--+  run_test


function build_tests()
{
   local make_args="MFEM_DIR=${MFEM_DIR} BLD=${test_exe_dir}"
   quoted_echo make -j $num_proc_build $make_args
   [[ -n "$dry_run" ]] || make -j $num_proc_build $make_args
}


function split3_power2()
{
   local outvar="$1" number="$2" n1="$3" n2="$4" n3="$5" nn=
   (( nn = 2**(number/3) ))
   (( n1 = n1*nn*(2**(number%3>0?1:0)) ))
   (( n2 = n2*nn*(2**(number%3>1?1:0)) ))
   (( n3 = n3*nn ))
   eval "$outvar=($n1 $n2 $n3)"
}


function make_mesh_file()
{
   [[ ${#mesh_nxyz[@]} -eq 3 ]] || {
      echo "Invalid mesh_nxyz variable."; exit 1; }
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
}


function run_test()
{
   set_mpi_options
   local common_args="-no-vis $mesh_opt -rs $ser_ref -rp $par_ref -o $sol_p"
   local all_args=()
   all_args=($common_args "${test_extra_args[@]}")
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

test_name=bp_main
# problem: 0 - diffusion, 1 - mass
problem=${problem:-0}
# quadrature type: 0 - Gauss, 1 - Gauss-Lobatto
ir_type=${ir_type:-0}
vdim=${vdim:-1}
vec_layout=${vec_layout:-}
# test id:     1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
sol_p_list=(   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)
ir_order_list=(0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0)
enabled_tests_def="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16"
# 'enabled_tests' can be set on the 'go.sh' command line
enabled_tests="${enabled_tests:-$enabled_tests_def}"
ser_ref=0
mesh_s_reduction_base=0     # used in run_tests()
mesh_s_reduction_limit=30   # used in run_tests()
mesh_max_elem=$(( 2**29 )) # 2^30 causes overflow in number of faces/edges
max_approx_dofs=$((5*(2**20)*num_proc_run))
# 'max_dofs' can be set on the 'go.sh' command line
max_approx_dofs="${max_dofs:-$max_approx_dofs}"
echo "max_approx_dofs: $max_approx_dofs"
# 'base_nxyz' can be set on the 'go.sh' command line
base_nxyz=(${base_nxyz[@]:-1 1 1})
base_n=$(( base_nxyz[0]*base_nxyz[1]*base_nxyz[2] ))

local n=$(( num_proc_run/base_n ))
proc_t=0
while (( n >= 2 )); do
   ((proc_t=proc_t+1))
   ((n=n/2))
done
(( num_proc_run == base_n*(2**proc_t) )) || {
   echo "Invalid number of processors: $num_proc_run"
   echo "The number of processors must be $base_n times a power of 2."
   return 1
}
split3_power2 proc_nxyz $proc_t "${base_nxyz[@]}"
test_extra_args=(-c "${proc_nxyz[*]}")

return 0

}


# Called by 'run_tests'
function set_test_params()
{
   [[ -n "$1" ]] && {
      local tid=$(( $1 - 1 ))
      sol_p="${sol_p_list[$tid]}"
      ir_order="${ir_order_list[$tid]}"
      local max_elems=$(( max_approx_dofs/((sol_p**3)*vdim) ))
      mesh_s="0"
      while (( base_n*(2**(mesh_s+1)) <= max_elems )); do
         (( mesh_s = mesh_s + 1 ))
      done
      (( mesh_s = mesh_s + mesh_s_shift ))
      par_ref="0"
      while (( mesh_s >= 3 )); do
         (( mesh_s = mesh_s - 3 ))
         (( par_ref = par_ref + 1 ))
      done
   }
   # echo "proc_t = $proc_t, mesh_s = $mesh_s, ser_ref = $ser_ref"
   if (( mesh_s + 3*(ser_ref+par_ref) < 0 )); then
      # invalid number of elements
      return 1
   fi
   local num_elem=$(( base_n*(2**(mesh_s+3*(ser_ref+par_ref))) ))
   if (( num_elem > mesh_max_elem )); then
      # too many elements
      return 1
   fi
   local approx_dofs=$(( num_elem*(sol_p**3)*vdim ))
   if (( approx_dofs > max_approx_dofs )); then
      # too many dofs
      return 1
   fi
   # adjust mesh_s as needed:
   if (( proc_t > mesh_s + 3*ser_ref )); then
      if (( proc_t <= mesh_s + 3*(ser_ref+par_ref) )); then
         while (( proc_t > mesh_s + 3*ser_ref )); do
            (( mesh_s = mesh_s + 3 ))
            (( par_ref = par_ref - 1 ))
         done
      else
         # number of processors is larger than number of elements
         return 1
      fi
   fi
   split3_power2 mesh_nxyz $mesh_s "${base_nxyz[@]}"
   make_mesh_file
   mesh_opt="-m $mesh_file"
}


# Called one time by the main function, 'build_and_run_tests'
function run_tests()
{
   local test_id= sft=
   for test_id; do
      for (( sft = mesh_s_reduction_base;
             sft <= mesh_s_reduction_limit; sft++ )); do
         (( mesh_s_shift = -sft ))
         set_test_params $test_id && run_test || {
            if (( mesh_s+3*(ser_ref+par_ref) >= 0 )); then
               printf "Skipped test #$test_id with number of mesh elements"
               echo " = $(( base_n*(2**(mesh_s+3*(ser_ref+par_ref))) ))."
            fi
         }
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


test_required_packages="metis hypre occa-fork mfem-engines"
