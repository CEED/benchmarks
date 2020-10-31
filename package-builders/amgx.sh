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

# Clone and build AMGX.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="amgx"
AMGX_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/amgx"
AMGX_DIR="$pkg_bld_dir"
pkg_var_prefix="amgx_"
amgx_branch="${amgx_branch:-main}"
pkg="AMGX"


function amgx_clone()
{
   pkg_repo_list=("git@github.com:NVIDIA/AMGX.git"
                  "https://github.com/NVIDIA/AMGX.git")
   pkg_git_branch="$amgx_branch"
   cd "$pkg_sources_dir" || return 1
   if [[ -d "$pkg_src_dir" ]]; then
      update_git_package
      return
   fi
   for pkg_repo in "${pkg_repo_list[@]}"; do
      echo "Cloning $pkg from $pkg_repo ..."
      git clone --recursive "$pkg_repo" "$pkg_src_dir" && return 0
   done
   echo "Could not successfully clone $pkg. Stop."
   return 1
}


function amgx_build()
{
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   elif [[ ! -d "$pkg_bld_dir" ]]; then
     # Start the build process from scratch.
     echo "Starting $pkg_src_dir build process from scratch" &&
     rm -rf "$pkg_bld_dir"
     mkdir -p "$pkg_bld_dir" || {
      echo "Error creating directory $pkg_bld_dir. Stop."
      return 1
     }
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      cd "$AMGX_SOURCE_DIR" && \
      cd "$pkg_bld_dir" && \
      cmake "$AMGX_SOURCE_DIR" \
         -DCMAKE_C_COMPILER="$MPICC" \
         -DCMAKE_C_FLAGS="$CFLAGS" \
         -DCMAKE_CXX_COMPILER="$MPICXX" \
         -DCMAKE_CXX_FLAGS="$CFLAGS -std=c++11" \
         -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_INSTALL_PREFIX="$AMGX_DIR" \
         -DCMAKE_VERBOSE_MAKEFILE=1 && \
      make -j $num_proc_build && \
      make install
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build successful."
   : > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   amgx_clone && get_package_git_version && amgx_build
}
