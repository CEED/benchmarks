#!/bin/bash

function grep_data()
{
  cd sin

  for i in `seq $min_order 1 $max_order`
  do
    cd lx$i
    for j in `seq $min_elem 1 $max_elem`
    do
      grep "case vec" lx$i/b$j/logfile >> sin.vec
      grep "case sca" lx$i/b$j/logfile >> sin.sca
    done
    cd ..
  done

  
  cd ..
}

function postprocess()
{
  echo 'Postprocessing ...'
  cd $OUT_DIR
  echo `pwd`
}
