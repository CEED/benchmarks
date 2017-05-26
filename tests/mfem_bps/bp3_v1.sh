# This file is part of CEED. For more details, see exascaleproject.org.


if [[ -z "$root_dir" ]]; then
   echo "This script ($0) should not be called directly. Stop."
   return 1
fi

# problem: 0 - diffusion, 1 - mass
problem=0
source ${root_dir}/tests/mfem_bps/bp1_v1.sh
