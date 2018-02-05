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

# Clone and build the CEED benchmarks from the Virginia Tech CEED team:
#    https://github.com/kswirydo/CEED-Ax

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="ceed-ax"
CEED_AX_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/ceed-ax"
CEED_AX_DIR="$pkg_bld_dir"
pkg_var_prefix="ceed_ax_"
pkg="CEED-Ax BPs"


function ceed_ax_clone()
{
   pkg_repo_list=("git@github.com:kswirydo/CEED-Ax.git"
                  "https://github.com/kswirydo/CEED-Ax.git")
   pkg_git_branch="master"
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


function ceed_ax_build_aux()
{
   local bp_dir=
   for bp_dir; do
      cd "$pkg_bld_dir/$bp_dir" && \
      make \
         OCCA_DIR="${OCCA_DIR}" \
         CXX="$MPICXX" \
         CXXFLAGS="$CFLAGS" \
         $OCCA_EXTRA_CONFIG \
         -j $num_proc_build || return 1
   done
   return 0
}


function ceed_ax_build()
{
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   elif [[ ! -d "$pkg_bld_dir" ]]; then
      echo "Cloning $pkg from $CEED_AX_SOURCE_DIR to OUT_DIR ..."
      cd "$OUT_DIR" && git clone "$CEED_AX_SOURCE_DIR" || {
         echo "Cloning $CEED_AX_SOURCE_DIR to OUT_DIR failed. Stop."
         return 1
      }
   fi
   if [[ -z "$OCCA_DIR" ]]; then
      echo "The required variable 'OCCA_DIR' is not set. Stop."
      return 1
   fi
   if [[ "${OCCA_BRANCH}" != "master" ]]; then
      echo "The package $pkg requires 'occa_branch=master'. Stop."
      return 1
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      ceed_ax_build_aux "BP10" "BP30" "BP35"
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build successful."
   print_variables "$pkg_var_prefix" \
      OCCA_DIR \
      > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   ceed_ax_clone && get_package_git_version && ceed_ax_build
}
