#!/bin/bash

function build_and_run_tests()
{
   export BP_ROOT="$root_dir"/tests/nek5000_bps
   export BENCH_ROOT="$root_dir"

   # Generate the boxes
   cd $BP_ROOT/boxes
   ./boxes.sh
   cd "$test_exe_dir"

   # Setup the sin version of the bp
   mkdir sin
   cd sin

   for i in `seq 2 1 4`
   do

      mkdir lx$i
      cp -r $BP_ROOT/boxes/b?? $BP_ROOT/SIZE $BP_ROOT/bp1/zsin.usr lx$i/

      # Set lx1 in SIZE file
      sed -i "s/lx1=[0-9]*/lx1=${i}/" lx$i/SIZE

      # Make the executable
      cd lx$i
      $BP_ROOT/makenek zsin
      cd ..

   done

   return 0
}
