# This file is part of CEED. For more details, see exascaleproject.org.

# Clone Nek5000.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
pkg_src_dir="Nek5000"
NEK5K_SOURCE_DIR="$pkg_sources_dir/$pkg_src_dir"
pkg_bld_dir="$OUT_DIR/nek5000"
MFEM_DIR="$pkg_bld_dir"
pkg="NEK5K"


function nek5k_clone()
{
   pkg_repo_list=("git@github.com:Nek5000/Nek5000.git"
                  "https://github.com/Nek5000/Nek5000.git")
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

function nek5k_build()
{
## Just build the requited tools: genbox and genmap
   cd $NEK5K_SOURCE_DIR

   cp bin/makenek ../../tests/nek5000_bps/ 

   cd tools
   ./maketools genmap
   ./maketools genbox

   return 0
}

function build_package()
{
   nek5k_clone && nek5k_build
}
