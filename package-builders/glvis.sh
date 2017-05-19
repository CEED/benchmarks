# This file is part of CEED. For more details, see exascaleproject.org.

# Clone and build GLVis.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   exit 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   exit 1
fi
pkg_src_dir="glvis"
GLVIS_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/glvis"
GLVIS_DIR="$pkg_bld_dir"


function glvis_clone()
{
   local pkg="GLVis"
   pkg_repo_list=("git@github.com:glvis/glvis.git"
                  "https://github.com/glvis/glvis.git")
   cd "$pkg_sources_dir" || return 1
   [[ -d "$pkg_src_dir" ]] && return 0
   for pkg_repo in "${pkg_repo_list[@]}"; do
      echo "Cloning $pkg from $pkg_repo ..."
      git clone "$pkg_repo" "$pkg_src_dir" && return 0
   done
   echo "Could not successfully clone $pkg. Stop."
   exit 1
}


function glvis_build()
{
   local pkg="GLVis"
   if [[ ! -d "$pkg_bld_dir" ]]; then
      cd "$OUT_DIR" && git clone "$GLVIS_SOURCE_DIR" || {
         echo "Cloning $GLVIS_SOURCE_DIR to OUT_DIR failed. Stop."
         exit 1
      }
   elif [[ -e "${pkg_bld_dir}_build_successful" ]]; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   fi
   if [[ -z "$MFEM_SERIAL_DIR" ]]; then
      echo "The required variable 'MFEM_SERIAL_DIR' is not set. Stop."
      exit 1
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      cd "$pkg_bld_dir" && \
      make \
         MFEM_DIR="$MFEM_SERIAL_DIR" \
         $GLVIS_EXTRA_CONFIG \
         -j $num_proc_build
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      exit 1
   }
   echo "Build succesful."
   : > "${pkg_bld_dir}_build_successful"
}


glvis_clone && glvis_build
