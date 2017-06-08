#!/bin/bash

function grep_data()
{
  cd sin

  for i in `seq $min_order 1 $max_order`
  do
    cd lx$i
    rm sin.vec sin.sca

    for j in `seq $min_elem 1 $max_elem`
    do
      grep "case vec" b$j/logfile >> sin.vec
      grep "case sca" b$j/logfile >> sin.sca
    done

    cd ..
  done

  cd ..
}

function plot_data()
{
  cd sin 

  gnuplot -e "start=$min_order; end=$max_order" \
     $root_dir"/tests"$test_up_dir/plot.gp || exit 1
}

function postprocess()
{
  # Get the min, max order and min, max elements
  . $root_dir"/tests"$test_up_dir/bp1.sh || exit 1
  configure_tests

  echo 'Postprocessing ...'
  cd $test_exe_dir
  grep_data

  # Plot with gnuplot
  plot_data
}