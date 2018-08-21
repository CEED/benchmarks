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

# Clone and build libCEED.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="libceed"
LIBCEED_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/libceed"
LIBCEED_DIR="$pkg_bld_dir"
libceed_branch="${libceed_branch:-master}"
LIBCEED_BRANCH="${libceed_branch}"
pkg_var_prefix="libceed_"
pkg="libCEED (branch $libceed_branch)"


function libceed_clone()
{
   pkg_repo_list=("git@github.com:CEED/libCEED.git"
                  "https://github.com/CEED/libCEED.git")
   pkg_git_branch="$libceed_branch"
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


function libceed_build()
{
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   fi
   # Start the build from scratch
   [[ -d "${pkg_bld_dir}" ]] && remove_package
   printf "Cloning libCEED from the package source directory "
   echo "$LIBCEED_SOURCE_DIR to the build directory inside OUT_DIR ..."
   cd "$OUT_DIR" && git clone "$LIBCEED_SOURCE_DIR" || {
      echo "Cloning failed. Stop."
      return 1
   }
   # Check for optional packages used by backends
   local OCCA_MAKE_OPTS=()
   if [[ -n "$OCCA_DIR" ]]; then
      OCCA_MAKE_OPTS=("OCCA_DIR=${OCCA_DIR}")
   else
      echo "${magenta}Warning: Building $pkg without OCCA ...${none}"
   fi
   # TODO: other packages: magma, ...
   libceed_verbose=${libceed_verbose:-1}
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      cd "$pkg_bld_dir" && \
      make -j $num_proc_build \
         V="$libceed_verbose" \
         CC="$MPICC" \
         FC="$MPIF77" \
         OPT="$CFLAGS" \
         "${OCCA_MAKE_OPTS[@]}"
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build successful."
   print_variables "$pkg_var_prefix" \
      LIBCEED_BRANCH OCCA_DIR \
      > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   libceed_clone && get_package_git_version && libceed_build
}