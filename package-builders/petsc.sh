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

# Clone and build PETSc.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="petsc"
PETSC_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/petsc"
PETSC_DIR="$pkg_bld_dir"
PETSC_ARCH="$short_config"
pkg_var_prefix="petsc_"
pkg="PETSc"


function petsc_clone()
{
   pkg_repo_list=("https://bitbucket.org/petsc/petsc.git"
                  "git@bitbucket.org:petsc/petsc.git")
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


function petsc_build()
{
   if package_build_is_good; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   elif [[ ! -d "$pkg_bld_dir" ]]; then
      echo "Cloning PETSc from $PETSC_SOURCE_DIR to OUT_DIR ..."
      cd "$OUT_DIR" && git clone "$PETSC_SOURCE_DIR" || {
         echo "Cloning $PETSC_SOURCE_DIR to OUT_DIR failed. Stop."
         return 1
      }
   fi
   if [[ -z "$HYPRE_DIR" ]]; then
      echo "The required variable 'HYPRE_DIR' is not set. Stop."
      return 1
   fi
   local PETSC_LAPACK_CONF=()
   if [[ -n "$LAPACK_LIB" ]]; then
      PETSC_LAPACK_CONF=("--with-blas-lib=$LAPACK_LIB"
                         "--with-lapack-lib=")
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      cd "$pkg_bld_dir" && \
      if [[ -e configure.log ]]; then
         git clean -df
         rm -f configure.log
      fi && \
      ./configure \
         --with-petsc-arch="$PETSC_ARCH" \
         --with-shared-libraries="0" \
         --with-cc="$MPICC" \
         --CFLAGS="$CFLAGS" \
         --with-cxx="$MPICXX" \
         --CXXFLAGS="$CFLAGS" \
         --with-fc="$MPIF77" \
         --FFLAGS="$FFLAGS" \
         --with-debugging="0" \
         --with-mpi="1" \
         "${PETSC_LAPACK_CONF[@]}" \
         --with-hypre="1" \
         --with-hypre-dir="$HYPRE_DIR/src/hypre" && \
      make -j $num_proc_build
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build successful."
   print_variables "$pkg_var_prefix" \
      HYPRE_DIR > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   petsc_clone && get_package_git_version && petsc_build
}
