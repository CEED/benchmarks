# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
# reserved. See files LICENSE and NOTICE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research is supported by the Exascale Computing Project
# (17-SC-20-SC), a collaborative effort of two U.S. Department of Energy
# organizations (Office of Science and the National Nuclear Security
# Administration) responsible for the planning and preparation of a capable
# exascale ecosystem, including software, applications, hardware, advanced
# system engineering and early testbed platforms, in support of the nation's
# exascale computing imperative.


if [[ -z "$bp_test" ]]; then
   echo "This script (bp-main.sh) should not be called directly."
   echo "Use one of the bp*.sh scripts instead. Stop."
   return 1
fi


function build_and_run_tests()
{
   local bp_exe="$DEALII_CEED_BPS_DIR/$bp_test"

   if [[ ! -x "$bp_exe" ]]; then
      echo "Invalid test: $bp_exe. Stop."
      return 1
   fi

   set_mpi_options

   $dry_run $mpi_run $bp_exe
}


test_required_packages="p4est dealii dealii-ceed-bps"
