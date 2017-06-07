# This file is part of CEED. For more details, see exascaleproject.org.

# Clone MFEM and build the parallel version using the occa-dev branch.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="mfem-occa"
MFEM_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/mfem-occa"
MFEM_DIR="$pkg_bld_dir"
pkg_version="$(git --git-dir=$MFEM_SOURCE_DIR/.git describe --long --abbrev=10 --tags)"
pkg="MFEM (occa-dev)"


function mfem_occa_clone()
{
   pkg_repo_list=("git@github.com:mfem/mfem.git"
                  "https://github.com/mfem/mfem.git")
   pkg_git_branch="occa-dev"
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


function mfem_occa_build()
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
   local METIS_5="NO"
   [[ "$METIS_VERSION" = "5" ]] && METIS_5="YES"
   if [[ -z "$OCCA_DIR" ]]; then
      echo "The required variable 'OCCA_DIR' is not set. Stop."
      return 1
   fi
   local SUNDIALS_MAKE_OPTS=()
   if [[ -n "$SUNDIALS_DIR" ]]; then
      SUNDIALS_MAKE_OPTS=(
         "MFEM_USE_SUNDIALS=YES"
         "SUNDIALS_DIR=$SUNDIALS_DIR")
   else
      echo "${magenta}Warning: Building $pkg without SUNDIALS ...${none}"
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      local num_nodes=1  # for 'make check' or 'make test'
      set_mpi_options    # for 'make check' or 'make test'
      cd "$pkg_bld_dir" && \
      make config \
         -f "$MFEM_SOURCE_DIR/makefile" \
         MFEM_USE_MPI=YES \
         MFEM_USE_OCCA=YES \
         $MFEM_EXTRA_CONFIG \
         MPICXX="$MPICXX" \
         CXXFLAGS="$CFLAGS" \
         HYPRE_DIR="$HYPRE_DIR/src/hypre" \
         METIS_DIR="$METIS_DIR" \
         MFEM_USE_METIS_5="$METIS_5" \
         OCCA_DIR="$OCCA_DIR" \
         "${SUNDIALS_MAKE_OPTS[@]}" \
         LDFLAGS="${LDFLAGS[*]}" \
         MFEM_MPIEXEC="${MPIEXEC:-mpirun}" \
         MFEM_MPIEXEC_NP="${MPIEXEC_OPTS} ${MPIEXEC_NP:--np}" && \
      make -j $num_proc_build
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build succesful."
   : > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   mfem_occa_clone && mfem_occa_build
}
