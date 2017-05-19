# This file is part of CEED. For more details, see exascaleproject.org.

# Clone MFEM and build the serial version without any dependencies.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   exit 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   exit 1
fi
pkg_src_dir="mfem"
MFEM_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/mfem-serial"
MFEM_SERIAL_DIR="$pkg_bld_dir"


function mfem_clone()
{
   local pkg="MFEM"
   pkg_repo_list=("git@github.com:mfem/mfem.git"
                  "https://github.com/mfem/mfem.git")
   cd "$pkg_sources_dir" || return 1
   [[ -d "$pkg_src_dir" ]] && return 0
   for pkg_repo in "${pkg_repo_list[@]}"; do
      echo "Cloning $pkg from $pkg_repo ..."
      git clone "$pkg_repo" "$pkg_src_dir" && return 0
   done
   echo "Could not successfully clone $pkg. Stop."
   exit 1
}


function mfem_serial_build()
{
   local pkg="mfem-serial"
   if [[ ! -d "$pkg_bld_dir" ]]; then
      mkdir -p "$pkg_bld_dir"
   elif [[ -e "${pkg_bld_dir}_build_successful" ]]; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      cd "$pkg_bld_dir" && \
      make config \
         -f "$MFEM_SOURCE_DIR/makefile" \
         $MFEM_EXTRA_CONFIG \
         CXX="$mpi_cxx" \
         CXXFLAGS="$CFLAGS" && \
      make -j $num_proc_build
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      exit 1
   }
   echo "Build succesful."
   : > "${pkg_bld_dir}_build_successful"
}


mfem_clone && mfem_serial_build
