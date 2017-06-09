#!/bin/bash

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

  cp $root_dir"/tests"$test_up_dir"/plot.gp" . 

  printf "set output \"%s_%s.png\"\n" "$1" "$2" >> plot.gp
  if [[ "$2" == "vec" ]]; then
    printf "set title \"Nek5000 - BP1 - %s - vector\" font \",16\"\n" "$1" >> plot.gp
  else
    printf "set title \"Nek5000 - BP1 - %s - scalar\" font \",16\"\n" "$1" >> plot.gp
  fi 

  if [[ "$min_order" != "$max_order" ]]; then
    printf "plot \"lx%d/%s.%s\" using 7:11 title 'lx%d' with linespoints,\\" "$min_order" "$1" "$2"  "$min_order"  >> plot.gp
    printf "\n" >> plot.gp

    start=$(( $min_order + 1 ))
    endt=$(( $max_order - 1 ))
    for i in `seq $start 1 $end`
    do
      printf "     \"lx%d/%s.%s\" using 7:11 title 'lx%d' with linespoints,\\" "$i" "$1" "$2" "$i"  >> plot.gp
      printf "\n" >> plot.gp
    done

    printf "     \"lx%d/%s.%s\" using 7:11 title 'lx%d' with linespoints" "$max_order" "$1" "$2" "$max_order"  >> plot.gp
  else
    printf "plot \"lx%d/%s.%s\" using 7:11 title 'lx%d' with linespoints" "$min_order" "$1" "$2"  "$min_order"  >> plot.gp
  fi

  gnuplot plot.gp

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
  grep_data zw 
 
  echo 'Plotting with gnuplot ...'
  plot_data zsin vec 
  plot_data zsin sca
  plot_data zw vec
  plot_data zw sca
}
