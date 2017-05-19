# This file is part of CEED. For more details, see exascaleproject.org.

# Download and build METIS v5.

if [[ -z "$root_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   exit 1
fi
METIS_VERSION="5"
. "$root_dir/package-builders/metis.sh"
