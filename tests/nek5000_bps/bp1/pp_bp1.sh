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
# testbed platforms, in support of the nation’s exascale computing imperative.

function grep_data()
{
  cd $1

  for i in `seq $min_order 1 $max_order`
  do
    cd lx$i
    rm $1.vec $1.sca

    for j in `seq $min_elem 1 $max_elem`
    do
      grep "case vec" b$j/logfile >> $1.vec
      grep "case sca" b$j/logfile >> $1.sca
    done

    cd ..
  done

  cd ..
}

function plot_data()
{
  cd $1

  rm plot_$2.gp

  printf "set xlabel \"DOFS\" font \",14\"\n" >> plot_$2.gp
  printf "set ylabel \"DOFS/s\" font \",14\"\n" >> plot_$2.gp
  printf "set logscale y\n" >> plot_$2.gp
  printf "set logscale x\n" >> plot_$2.gp
  printf "set term png\n\n" >> plot_$2.gp

  printf "set output \"%s_%s.png\"\n" "$1" "$2" >> plot_$2.gp
  if [[ "$2" == "vec" ]]; then
    printf "set title \"Nek5000 - BP1 - %s - vector\" font \",16\"\n" "$1" >> plot_$2.gp
  else
    printf "set title \"Nek5000 - BP1 - %s - scalar\" font \",16\"\n" "$1" >> plot_$2.gp
  fi

  if [[ "$min_order" != "$max_order" ]]; then
    printf "plot \"lx%d/%s.%s\" using 7:11 title 'lx%d' with linespoints,\\" "$min_order" "$1" "$2"  "$min_order"  >> plot_$2.gp
    printf "\n" >> plot_$2.gp

    start=$(( $min_order + 1 ))
    end=$(( $max_order - 1 ))
    for i in `seq $start 1 $end`
    do
      printf "     \"lx%d/%s.%s\" using 7:11 title 'lx%d' with linespoints,\\" "$i" "$1" "$2" "$i"  >> plot_$2.gp
      printf "\n" >> plot_$2.gp
    done

    printf "     \"lx%d/%s.%s\" using 7:11 title 'lx%d' with linespoints" "$max_order" "$1" "$2" "$max_order"  >> plot_$2.gp
  else
    printf "plot \"lx%d/%s.%s\" using 7:11 title 'lx%d' with linespoints" "$min_order" "$1" "$2"  "$min_order"  >> plot_$2.gp
  fi

  gnuplot plot_$2.gp

  cd ..
}

function postprocess()
{
  # Get the min, max order and min, max elements
  . $root_dir"/tests"$test_up_dir/bp1.sh || exit 1
  configure_tests

  echo 'Postprocessing ...'
  cd $test_exe_dir
  grep_data zsin
  # w tests are not run by default
#  grep_data zw

  echo 'Plotting with gnuplot ...'
  plot_data zsin vec
  plot_data zsin sca
  # w tests are not run by default
#  plot_data zw vec
#  plot_data zw sca
}
