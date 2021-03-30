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

# Clone and build RAJA.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="raja"
RAJA_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/raja"
RAJA_DIR="$pkg_bld_dir/install"
pkg_var_prefix="raja_"
pkg="RAJA"


function raja_clone()
{
   pkg_repo_list=("git@github.com:LLNL/raja.git"
                  "https://github.com/LLNL/raja.git")
   pkg_git_branch="main"
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


function raja_build()
{
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   fi
   # Start the build process from scratch.
   rm -rf "$pkg_bld_dir"
   local raja_cuda_opts=()
   if [[ -n "$CUDA_ENABLED" ]]; then
      raja_cuda_opts=(
         -DENABLE_CUDA="ON"
         -DCUDA_TOOLKIT_ROOT_DIR="$cuda_home"
         -DCMAKE_CUDA_COMPILER="$cuda_home/bin/nvcc"
         -DCUDA_ARCH="${cuda_arch:-sm_35}")
   else
      raja_cuda_opts=(-DENABLE_CUDA="OFF")
      echo "${magenta}INFO: Building $pkg without CUDA ...${none}"
   fi
   local raja_openmp="OFF"
   if [[ -n "$OMP_ENABLED" ]]; then
      raja_openmp="ON"
   else
      echo "${magenta}INFO: Building $pkg without OpenMP ...${none}"
   fi
   local raja_hip_opts=()
   if [[ -n "$HIP_ENABLED" ]]; then
      raja_hip_opts=(
        -DENABLE_HIP="YES"
        -DHIP_ROOT_DIR="$hip_home"
        -DHIP_RUNTIME_INCLUDE_DIRS="$hip_home/include"
        -DHIP_HIPCC_FLAGS="--amdgpu-target=$hip_arch"
      )
   else
      raja_hip_opts=(-DENABLE_HIP="OFF")
      echo "${magenta}INFO: Building $pkg without HIP ...${none}"
   fi
   mkdir -p "$pkg_bld_dir" || {
      echo "Error creating directory $pkg_bld_dir. Stop."
      return 1
   }
   # Note: RAJA does not seem to use the values of variables like
   #       CMAKE_CXX_FLAGS, CMAKE_CXX_FLAGS_RELEASE, etc, at least when CUDA
   #       support is enabled.
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      cd "$RAJA_SOURCE_DIR" && \
      git submodule update && \
      cd "$pkg_bld_dir" && \
      cmake "$RAJA_SOURCE_DIR" \
         -DCMAKE_C_COMPILER="$MPICC" \
         -DCMAKE_C_FLAGS="$CFLAGS" \
         -DCMAKE_CXX_COMPILER="$MPICXX" \
         -DCMAKE_CXX_FLAGS="$CFLAGS" \
         -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_INSTALL_PREFIX="$RAJA_DIR" \
         -DENABLE_TESTS=OFF \
         -DENABLE_EXAMPLES=OFF \
         -DENABLE_EXERCISES=OFF \
         -DENABLE_MPI=ON \
         "${raja_cuda_opts[@]}" \
         "${raja_hip_opts[@]}" \
         -DENABLE_OPENMP="$raja_openmp" \
         -DCMAKE_VERBOSE_MAKEFILE=1 && \
      make -j $num_proc_build && \
      make install
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build successful."
   print_variables "$pkg_var_prefix" \
      CUDA_ENABLED cuda_home OMP_ENABLED omp_flag \
      > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   raja_clone && get_package_git_version && raja_build
}
