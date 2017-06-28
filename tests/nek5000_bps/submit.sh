#!/bin/bash

# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-XXXXXX. All Rights
# reserved. See file LICENSE for details.
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

#SBATCH -o out.file
#SBATCH -e error.file

echo $1        >  SESSION.NAME
echo `pwd`'/' >>  SESSION.NAME
touch $1.rea
rm -f logfile
rm -f ioinfo
mv $1.log.$2 $1.log1.$2 2>/dev/null
mv $1.sch $1.sch1       2>/dev/null
# echo "Executing: $mpi_run ./nek5000 > $1.log.$2"
echo "Executing: $mpi_run ./nek5000"
echo "In directory: $PWD"
# $mpi_run ./nek5000 > $1.log.$2
$mpi_run ./nek5000
sleep 2
# ln $1.log.$2 logfile
