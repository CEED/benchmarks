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

# Clone and build the development version of deal.II.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="dealii"
DEALII_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/dealii"
DEALII_DIR="$pkg_bld_dir"
pkg_var_prefix="dealii_"
pkg="deal.II"


function dealii_clone()
{
   pkg_repo_list=("git@github.com:dealii/dealii.git"
                  "https://github.com/dealii/dealii.git")
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


function dealii_build()
{
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   elif [[ ! -d "$pkg_bld_dir" ]]; then
      mkdir -p "$pkg_bld_dir"
   fi
   if [[ -z "$P4EST_DIR" ]]; then
      echo "The required variable 'P4EST_DIR' is not set. Stop."
      return 1
   fi
   local LAPACK_CMAKE_OPTS=("-DDEAL_II_WITH_LAPACK=OFF")
   if [[ -n "$LAPACK_LIB" ]]; then
      LAPACK_CMAKE_OPTS=(
         "-DDEAL_II_WITH_LAPACK=ON"
         "-DLAPACK_LIBRARIES=$LAPACK_LIB")
      array_union LDFLAGS "$LAPACK_LIB"
   else
      echo "${magenta}Warning: Building $pkg without LAPACK ...${none}"
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      mkdir -p "$pkg_bld_dir"/build && cd "$pkg_bld_dir"/build && \
      cmake "$DEALII_SOURCE_DIR" \
         -DCMAKE_BUILD_TYPE="Release" \
         -DCMAKE_VERBOSE_MAKEFILE=1 \
         -DCMAKE_C_COMPILER="$MPICC" \
         -DCMAKE_C_FLAGS_RELEASE="$CFLAGS" \
         -DCMAKE_CXX_COMPILER="$MPICXX" \
         -DCMAKE_CXX_FLAGS_RELEASE="$CFLAGS" \
         -DCMAKE_CXX_FLAGS="$NATIVE_CFLAG" \
         -DCMAKE_INSTALL_PREFIX="$pkg_bld_dir" \
         -DDEAL_II_WITH_MPI="ON" \
         "${LAPACK_CMAKE_OPTS[@]}" \
         -DDEAL_II_WITH_P4EST="ON" \
         -DP4EST_DIR="$P4EST_DIR" \
         -DDEAL_II_WITH_BOOST="OFF" \
         -DDEAL_II_WITH_THREADS="OFF" && \
      make -j $num_proc_build && \
      make install && \
      cd .. && rm -rf build
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build successful."
   print_variables "$pkg_var_prefix" \
      CFLAGS LAPACK_LIB P4EST_DIR > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   dealii_clone && get_package_git_version && dealii_build
}
