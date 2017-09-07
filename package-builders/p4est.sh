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

# Clone and build the development version of p4est.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="p4est"
P4EST_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/p4est"
P4EST_DIR="$pkg_bld_dir"
pkg_var_prefix="p4est_"
pkg="p4est"


function p4est_clone()
{
   pkg_repo_list=("git@github.com:cburstedde/p4est.git"
                  "https://github.com/cburstedde/p4est.git")
   pkg_git_branch="master"
   cd "$pkg_sources_dir" || return 1
   if [[ -d "$pkg_src_dir" ]]; then
      update_git_package && \
      cd "$pkg_src_dir" && git submodule init && git submodule update
      return
   fi
   for pkg_repo in "${pkg_repo_list[@]}"; do
      echo "Cloning $pkg from $pkg_repo ..."
      git clone "$pkg_repo" "$pkg_src_dir" && \
      cd "$pkg_src_dir" && git submodule init && git submodule update && \
      return 0
   done
   echo "Could not successfully clone $pkg. Stop."
   return 1
}


function p4est_build()
{
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   elif [[ ! -d "$pkg_bld_dir" ]]; then
      mkdir -p "$pkg_bld_dir"
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      cd "$P4EST_SOURCE_DIR" && \
      echo "--- bootstrap ----------------------------------------------" && \
      ./bootstrap && \
      echo "--- configure ----------------------------------------------" && \
      mkdir -p "$pkg_bld_dir"/build && cd "$pkg_bld_dir"/build && \
      $P4EST_SOURCE_DIR/configure \
         --enable-mpi \
         --enable-shared \
         --disable-vtk-binary \
         --without-blas \
         --prefix="$pkg_bld_dir" \
         CC="$MPICC" \
         CFLAGS="$CFLAGS" \
         CPPFLAGS="-DSC_LOG_PRIORITY=SC_LP_ESSENTIAL" && \
      echo "--- make ---------------------------------------------------" && \
      make -j $num_proc_build && \
      echo "--- make install -------------------------------------------" && \
      make install && \
      cd .. && rm -rf build
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build successful."
   print_variables "$pkg_var_prefix" \
      CFLAGS > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   p4est_clone && get_package_git_version && p4est_build
}
