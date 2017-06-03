# This file is part of CEED. For more details, see exascaleproject.org.

# Download and build METIS v4/v5.

if [[ -z "$pkg_sources_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi
if [[ -z "$OUT_DIR" ]]; then
   echo "The variable 'OUT_DIR' is not set. Stop."
   return 1
fi
METIS_VERSION="${METIS_VERSION:-4}"
pkg_src_url_base="http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis"
if [[ "$METIS_VERSION" = "4" ]]; then
   METIS_FULL_VERSION="4.0.3"
else
   METIS_FULL_VERSION="5.1.0"
fi
pkg_src_dir="metis-$METIS_FULL_VERSION"
pkg_src_url="$pkg_src_url_base/${pkg_src_dir}.tar.gz"
METIS_SOURCE_FILE="$pkg_sources_dir/${pkg_src_dir}.tar.gz"
pkg_bld_dir="$OUT_DIR/$pkg_src_dir"
METIS_DIR="$pkg_bld_dir"
pkg="METIS v$METIS_FULL_VERSION"


function metis_download()
{
   cd "$pkg_sources_dir" || return 1
   [[ -e "$METIS_SOURCE_FILE" ]] && return 0
   echo "Downloading $pkg from $pkg_src_url ..."
   if command -v wget &> /dev/null; then
      wget "$pkg_src_url" && return 0
   elif command -v curl &> /dev/null; then
      curl -O "$pkg_src_url" && return 0
   fi
   echo "Could not download $pkg. Stop."
   return 1
}


function metis_build()
{
   if [[ ! -d "$pkg_bld_dir" ]]; then
      cd "$OUT_DIR" && tar zxf "$METIS_SOURCE_FILE" || {
         echo "Error extracting \"$METIS_SOURCE_FILE\". Stop."
         return 1
      }
   elif [[ -e "${pkg_bld_dir}_build_successful" ]]; then
      echo "Using successfully built $pkg from OUT_DIR."
      return 0
   fi
   echo "Building $pkg, sending output to ${pkg_bld_dir}_build.log ..." && {
      if [[ "$METIS_VERSION" = "4" ]]; then
         cd "$pkg_bld_dir" && \
         make -C Lib realclean && \
         make -C Lib CC="$MPICC" OPTFLAGS="$CFLAGS" -j $num_proc_build
      else
         cd "$pkg_bld_dir/build" && \
         if [[ ! -e Makefile ]]; then \
            if ! command -v cmake &> /dev/null; then
               echo "Building $pkg requires cmake (not found). Stop."
               return 1
            fi
            cmake .. \
               -DCMAKE_VERBOSE_MAKEFILE=1 \
               -DGKLIB_PATH="$pkg_bld_dir/GKlib" \
               -DCMAKE_C_COMPILER="$MPICC" \
               -DCMAKE_CXX_COMPILER="$MPICXX" \
               -DSHARED="0" \
               -DBUILD_SHARED_LIBS="0" \
               -DMETIS_INSTALL="ON" \
               -DCMAKE_INSTALL_NAME_DIR="$pkg_bld_dir/lib" \
               -DCMAKE_INSTALL_RPATH="$pkg_bld_dir/lib" \
               -DCMAKE_INSTALL_PREFIX="$pkg_bld_dir"
         fi && \
         make metis -j $num_proc_build && \
         mkdir -p ../lib && cp -p libmetis/*metis* ../lib
      fi
   } &> "${pkg_bld_dir}_build.log" || {
      echo " ... building $pkg FAILED, see log for details."
      return 1
   }
   echo "Build successful."
   : > "${pkg_bld_dir}_build_successful"
}


function build_package()
{
   metis_download && metis_build
}
