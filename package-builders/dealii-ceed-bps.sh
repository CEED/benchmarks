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

# Clone and build the deal.II CEED benchmarks from:
#    https://github.com/kronbichler/ceed_benchmarks_dealii

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="dealii-ceed-bps"
DEALII_CEED_BPS_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/dealii-ceed-bps"
DEALII_CEED_BPS_DIR="$pkg_bld_dir"
pkg_var_prefix="dealii_ceed_bps_"
pkg="deal.II-CEED-BPs"


function dealii_ceed_bps_clone()
{
   pkg_repo_list=("git@github.com:kronbichler/ceed_benchmarks_dealii.git"
                  "https://github.com/kronbichler/ceed_benchmarks_dealii.git")
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


function dealii_ceed_bps_build_aux()
{
   local bp_dir=
   for bp_dir; do
      mkdir -p "$pkg_bld_dir/$bp_dir" && \
      cd "$pkg_bld_dir/$bp_dir" && \
      cmake "$DEALII_CEED_BPS_SOURCE_DIR/$bp_dir" \
         -DCMAKE_BUILD_TYPE="Release" \
         -DCMAKE_VERBOSE_MAKEFILE=1 \
         -DCMAKE_C_COMPILER="$MPICC" \
         -DCMAKE_C_FLAGS_RELEASE="$CFLAGS" \
         -DCMAKE_CXX_COMPILER="$MPICXX" \
         -DCMAKE_CXX_FLAGS_RELEASE="$CFLAGS" \
         -DCMAKE_CXX_FLAGS="$NATIVE_CFLAG" \
         -DDEAL_II_DIR="$DEALII_DIR" && \
      make -j $num_proc_build || return 1
   done
   return 0
}


function dealii_ceed_bps_build()
{
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   elif [[ ! -d "$pkg_bld_dir" ]]; then
      mkdir -p "$pkg_bld_dir"
   fi
   if [[ -z "$DEALII_DIR" ]]; then
      echo "The required variable 'DEALII_DIR' is not set. Stop."
      return 1
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      dealii_ceed_bps_build_aux "bp1" "bp2"
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build successful."
   : > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   dealii_ceed_bps_clone && get_package_git_version && dealii_ceed_bps_build
}
