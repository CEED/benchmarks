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

# Clone and build Laghos.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="laghos"
LAGHOS_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/laghos"
LAGHOS_DIR="$pkg_bld_dir"
# The laghos_branch can be set at the command line of the go.sh call
laghos_branch="${laghos_branch:-master}"
pkg_var_prefix="laghos_"
pkg="Laghos (branch ${laghos_branch})"


function laghos_clone()
{
   pkg_repo_list=("git@github.com:CEED/Laghos.git"
                  "https://github.com/CEED/Laghos.git")
   pkg_git_branch="${laghos_branch}"
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


function laghos_build()
{
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   elif [[ ! -d "$pkg_bld_dir" ]]; then
      echo "Cloning Laghos from $LAGHOS_SOURCE_DIR to OUT_DIR ..."
      cd "$OUT_DIR" && git clone "$LAGHOS_SOURCE_DIR" || {
         echo "Cloning $LAGHOS_SOURCE_DIR to OUT_DIR failed. Stop."
         return 1
      }
   fi
   if [[ -z "$MFEM_DIR" ]]; then
      echo "The required variable 'MFEM_DIR' is not set. Stop."
      return 1
   fi
   if [[ -z "$laghos_skip_main" ]]; then
      echo "Building $pkg (main version), sending output to" \
           " ${pkg_bld_dir}_build.log ..." && {
         cd "$pkg_bld_dir" && \
         make \
            MFEM_DIR="$MFEM_DIR" \
            CONFIG_MK="$MFEM_DIR/share/mfem/config.mk" \
            TEST_MK="$MFEM_DIR/share/mfem/test.mk" \
            -j $num_proc_build
      } &> "${pkg_bld_dir}_build.log" || {
         echo " ... building $pkg FAILED, see log for details."
         return 1
      }
   else
      echo "Building $pkg (main version) SKIPPED"
      : > "${pkg_bld_dir}_build.log"
   fi
   if [[ -n "$CUDA_ENABLED" ]]; then
      local mpi_exe="$(which $MPICXX)"
      local mpi_home="$(dirname $mpi_exe)/.."
      local nvcc_flags="-x=cu -std=c++11 -m64 --restrict -Xcompiler -Wall"
      nvcc_flags+=" -arch=${cuda_arch:-sm_60}"
      local nvcc_libs="-Wl,-rpath,$cuda_home/lib64 -L$cuda_home/lib64"
      nvcc_libs+=" -lcuda -lcudart -lcudadevrt -lnvToolsExt"
      echo "Building $pkg (CUDA version), appending output to" \
           "${pkg_bld_dir}_build.log ..." && {
         cd "$pkg_bld_dir/cuda" && \
         make \
            CXX="nvcc -ccbin $MPICXX" \
            CXXFLAGS="-Xcompiler=\"$CFLAGS\" -I$MFEM_DIR/include/mfem" \
            MFEM_DIR="$MFEM_DIR" \
            CONFIG_MK="$MFEM_DIR/share/mfem/config.mk" \
            TEST_MK="$MFEM_DIR/share/mfem/test.mk" \
            NVCC_CXXFLAGS="$nvcc_flags" \
            NVCC_LIBS="$nvcc_libs" \
            -j $num_proc_build
      } &>> "${pkg_bld_dir}_build.log" || {
         echo " ... building $pkg FAILED, see log for details."
         return 1
      }
   fi
   echo "Build successful."
   print_variables "$pkg_var_prefix" \
      laghos_branch \
      MFEM_DIR laghos_skip_main CUDA_ENABLED cuda_home \
      > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   laghos_clone && get_package_git_version && laghos_build
}
