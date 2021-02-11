# Copyright (c) 2021, Lawrence Livermore National Security, LLC. Produced at
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

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
if [[ -z "$MFEM_SOURCE_DIR" ]]; then
   echo "mfem must be built"
   return 1
fi
pkg_src_dir="solvers"
SOLVERS_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/solvers"
mfem_solvers_branch="${mfem_solvers_branch:-main}"
pkg_var_prefix="solvers_"
pkg="mfem solvers bp (branch $mfem_solvers_branch)"

function solvers_clone()
{
   pkg_repo_list=("git@github.com:ceed/solvers.git"
                  "https://github.com/ceed/solvers.git")
   pkg_git_branch="${mfem_solvers_branch:-master}"
   cd "$pkg_sources_dir" || return 1
   if [[ -d "$pkg_src_dir" ]]; then
      update_git_package
      return
   fi
   for pkg_repo in "${pkg_repo_list[@]}"; do
      echo "Cloning $pkg from $pkg_repo ..."
      git clone "$pkg_repo" "$pkg_src_dir" && return 0
   done
   echo "Could not successfully clone $pkg. Stop."
   return 1

}

function mfem_solvers_build()
{
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   elif [[ ! -d "$pkg_bld_dir" ]]; then
      echo "Cloning solvers from $SOLVERS_SOURCE_DIR to OUT_DIR ..."
      cd "$OUT_DIR" && git clone "$SOLVERS_SOURCE_DIR" || {
         echo "Cloning $SOLVERS_SOURCE_DIR to OUT_DIR failed. Stop."
         return 1
      }
   fi
   cd "$pkg_bld_dir/mfem-solver-bp"
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      cd "$pkg_bld_dir/mfem-solver-bp" && \
      (echo "MFEM_DIR=$OUT_DIR/mfem" > Make.user) && \
      if [[ "$MFEM_BRANCH" == "split-nonconforming-tet" ]]; then
         echo "USER_CXXFLAGS=-DMFEM_SIMPLEX_LOR" >> Make.user
      fi && \
      make -j $num_proc_build
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build successful."
}


function build_package()
{
   solvers_clone && get_package_git_version && mfem_solvers_build
}

