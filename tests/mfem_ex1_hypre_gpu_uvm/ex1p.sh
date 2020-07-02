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

problem=${problem:-0}

function build_tests()
{
   local make_opts=(
      "-j" "$num_proc_build"
      "MFEM_DIR=$MFEM_DIR"
      "BLD=$test_exe_dir/")
   quoted_echo make "${make_opts[@]}"
   [[ -n "$dry_run" ]] || make "${make_opts[@]}"
}


function run_tests()
{
   local test_name="ex1_hypre_gpu_uvm"
   set_mpi_options
   # 'min_p' can be set on the command line
   local l_min_p=${min_p:-1}
   # 'max_p' can be set on the command line
   local l_max_p=${max_p:-8}
   # 'max_dofs_proc' can be set on the command line
   local l_max_dofs_proc=${max_dofs_proc:-2100000}
   # 'mfem_devs' can be set on the command line
   local l_mfem_devs=(${mfem_devs:-cuda:uvm})
   local dim=3
   local args=
   for dev in ${l_mfem_devs[@]}; do
      for ((p = l_min_p; p <= l_max_p; p++)) do
         for ((l = 0; (p**dim)*(2**l) <= l_max_dofs_proc; l++)) do
            args=(-o $p -l $l -d $dev)
            if [ -z "$dry_run" ]; then
               echo "Running test:"
               quoted_echo $mpi_run ./$test_name "${args[@]}"
               $mpi_run ./$test_name "${args[@]}" || \
                  printf "\nError in the test, error code: $?\n\n"
            else
               $dry_run $mpi_run ./$test_name "${args[@]}"
            fi
         done
      done
   done
}


function build_and_run_tests()
{
   $dry_run cd "$test_dir"

   build_tests || return 1
   echo

   [[ -n "$build_only" ]] && return

   $dry_run cd "$test_exe_dir"
   run_tests
}


mfem_branch=${mfem_branch:-hypre-cuda-dev}
libceed_branch=${libceed_branch:-master}

# Uncomment the next line to enable 64-bit HYPRE_Int:
# hypre_big_int=1

packages=${packages:-cuda metis hypre-gpu-uvm mfem}

test_required_packages=${packages}
