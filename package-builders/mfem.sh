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

# Clone MFEM and build the parallel version.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="mfem"
MFEM_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/mfem"
MFEM_DIR="$pkg_bld_dir/install"
# 'mfem_branch' can be set at the command line of the go.sh call
mfem_branch="${mfem_branch:-master}"
MFEM_BRANCH="${mfem_branch}"
MFEM_DEBUG="${mfem_debug:+YES}"
pkg_var_prefix="mfem_"
pkg="MFEM (branch $mfem_branch)"


function mfem_clone()
{
   pkg_repo_list=("git@github.com:mfem/mfem.git"
                  "https://github.com/mfem/mfem.git")
   pkg_git_branch="${mfem_branch:-master}"
   cd "$pkg_sources_dir" || return 1
   if [[ -d "$pkg_src_dir" ]]; then
      update_git_package
      return
   fi
   for pkg_repo in "${pkg_repo_list[@]}"; do
      echo "Cloning $pkg from $pkg_repo ..."
      git clone "$pkg_repo" "$pkg_src_dir" && \
      cd "$pkg_src_dir" && \
      git checkout "$pkg_git_branch" && \
      return 0
   done
   echo "Could not successfully clone $pkg. Stop."
   return 1
}


function mfem_build()
{
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   elif [[ ! -d "$pkg_bld_dir" ]]; then
      mkdir -p "$pkg_bld_dir"
   fi
   if [[ -z "$HYPRE_DIR" ]]; then
      echo "The required variable 'HYPRE_DIR' is not set. Stop."
      return 1
   fi
   if [[ -z "$METIS_DIR" ]]; then
      echo "The required variable 'METIS_DIR' is not set. Stop."
      return 1
   fi
   local cxx11_flag="${CXX11FLAG:--std=c++11}"
   local optim_flags="$cxx11_flag $CFLAGS"
   local xcompiler=""
   local METIS_5="NO"
   [[ "$METIS_VERSION" = "5" ]] && METIS_5="YES"
   local SIMD_MAKE_OPTS=()
   if [[ -n "$SIMD_ENABLED" ]]; then
      SIMD_MAKE_OPTS=("MFEM_USE_SIMD=YES")
   else
      echo "${magenta}INFO: Building $pkg without SIMD ...${none}"
   fi
   local CUDA_MAKE_OPTS=()
   if [[ -n "$CUDA_ENABLED" ]]; then
      CUDA_MAKE_OPTS=(
         "MFEM_USE_CUDA=YES"
         "CUDA_CXX=$cuda_home/bin/nvcc"
         "CUDA_ARCH=${cuda_arch:-sm_70}")
      xcompiler="-Xcompiler="
      optim_flags="$cxx11_flag $xcompiler\"$CFLAGS\""
   else
      echo "${magenta}INFO: Building $pkg without CUDA ...${none}"
   fi
   local HIP_MAKE_OPTS=()
   if [[ -n "$HIP_ENABLED" ]]; then
      HIP_MAKE_OPTS=("MFEM_USE_HIP=YES"
                     "HIP_ARCH=${hip_arch}")
      echo "${cyan}INFO: Building $pkg with HIP ...${none}"
   else
      echo "${magenta}INFO: Building $pkg without HIP ...${none}"
   fi
   local OCCA_MAKE_OPTS=()
   if [[ -n "$OCCA_DIR" ]]; then
      OCCA_MAKE_OPTS=(
         "MFEM_USE_OCCA=YES"
         "OCCA_DIR=$OCCA_DIR")
   else
      echo "${magenta}INFO: Building $pkg without OCCA ...${none}"
   fi
   local RAJA_MAKE_OPTS=()
   if [[ -n "$RAJA_DIR" ]]; then
      RAJA_MAKE_OPTS=(
         "MFEM_USE_RAJA=YES"
         "RAJA_DIR=$RAJA_DIR")
   else
      echo "${magenta}INFO: Building $pkg without RAJA ...${none}"
   fi
   local AMGX_MAKE_OPTS=()
   if [[ -n "$AMGX_DIR" ]]; then
      AMGX_MAKE_OPTS=(
         "MFEM_USE_AMGX=YES"
         "AMGX_DIR=$AMGX_DIR")
   else
      echo "${magenta}INFO: Building $pkg without AMGX ...${none}"
   fi
   local OMP_MAKE_OPTS=()
   if [[ -n "$OMP_ENABLED" ]]; then
      OMP_MAKE_OPTS=(
         "MFEM_USE_OPENMP=YES"
         "OPENMP_OPT=$xcompiler\"$omp_flag\"")
   else
      echo "${magenta}INFO: Building $pkg without OpenMP ...${none}"
   fi
   local LIBCEED_MAKE_OPTS=()
   if [[ -n "$LIBCEED_DIR" ]]; then
      LIBCEED_MAKE_OPTS=(
         "MFEM_USE_CEED=YES"
         "CEED_DIR=$LIBCEED_DIR")
   else
      echo "${magenta}INFO: Building $pkg without libCEED ...${none}"
   fi
   local SUNDIALS_MAKE_OPTS=()
   if [[ -n "$SUNDIALS_DIR" ]]; then
      SUNDIALS_MAKE_OPTS=(
         "MFEM_USE_SUNDIALS=YES"
         "SUNDIALS_DIR=$SUNDIALS_DIR")
   else
      echo "${magenta}INFO: Building $pkg without SUNDIALS ...${none}"
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      local num_nodes=1  # for 'make check' or 'make test'
      set_mpi_options    # for 'make check' or 'make test'
      cd "$pkg_bld_dir" && \
      make config \
         -f "$MFEM_SOURCE_DIR/makefile" \
         PREFIX="$MFEM_DIR" \
         MFEM_USE_MPI=YES \
         ${mfem_debug:+MFEM_DEBUG=YES} \
         $MFEM_EXTRA_CONFIG \
         MPICXX="$MPICXX" \
         OPTIM_FLAGS="$optim_flags" \
         HYPRE_DIR="$HYPRE_DIR/src/hypre" \
         METIS_DIR="$METIS_DIR" \
         MFEM_USE_METIS_5="$METIS_5" \
         "${SIMD_MAKE_OPTS[@]}" \
         "${CUDA_MAKE_OPTS[@]}" \
         "${HIP_MAKE_OPTS[@]}" \
         "${OCCA_MAKE_OPTS[@]}" \
         "${RAJA_MAKE_OPTS[@]}" \
         "${AMGX_MAKE_OPTS[@]}" \
         "${OMP_MAKE_OPTS[@]}" \
         "${LIBCEED_MAKE_OPTS[@]}" \
         "${SUNDIALS_MAKE_OPTS[@]}" \
         LDFLAGS="${LDFLAGS[*]}" \
         MFEM_MPIEXEC="${MPIEXEC:-mpirun}" \
         MFEM_MPIEXEC_NP="${MPIEXEC_OPTS} ${MPIEXEC_NP:--np}" && \
      make info && \
      make -j $num_proc_build && \
      make install
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build successful."
   print_variables "$pkg_var_prefix" \
      MFEM_BRANCH MFEM_DEBUG \
      HYPRE_DIR METIS_DIR METIS_VERSION SIMD_ENABLED CUDA_ENABLED cuda_home \
      HIP_ENABLED hip_home OCCA_DIR RAJA_DIR OMP_ENABLED omp_flag \
      LIBCEED_DIR SUNDIALS_DIR \
      > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   mfem_clone && get_package_git_version && mfem_build
}
