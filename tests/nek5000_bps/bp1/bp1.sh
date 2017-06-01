#!/bin/bash

# Helper functions
function log2() 
{
    local x=0
    for (( y=$1-1 ; $y > 0; y >>= 1 )) ; do
        let x=$x+1
    done
    echo $x
}

function xyz() 
{
  prod=$1
  
  if [ $((prod%3)) -eq 0 ]
  then
    nez=$((prod/3)) 
    ney=$((prod/3)) 
    nex=$((prod/3)) 
  elif [ $((prod%3)) -eq 1 ] 
  then 
    nez=$((prod/3 +1)) 
    ney=$((prod/3)) 
    nex=$((prod/3)) 
  elif [ $((prod%3)) -eq 2 ] 
  then
    nez=$((prod/3 +1)) 
    ney=$((prod/3 +1)) 
    nex=$((prod/3)) 
  fi

  nex=$((2**nex))
  ney=$((2**ney))
  nez=$((2**nez))

  echo "$nex $ney $nez"
}

function genbb()
{
  cp $1.box ttt.box
  $NEK5K_DIR/bin/genbox <<EOF
ttt.box
EOF
$NEK5K_DIR/bin/genmap << EOF
box
.1
EOF
  $NEK5K_DIR/bin/mvn box $1

  rm ttt.box
  mv box.rea $1.rea

  return 0
}

function generate_boxes()
{
  # Run thorugh the box sizes
  for i in `seq $min_elem 1 $max_elem`
  do

    # Set the number of elements in box file.

    xyz=$(xyz $i)
    nex=$( echo $xyz | cut -f 1 -d ' ' )
    ney=$( echo $xyz | cut -f 2 -d ' ' )
    nez=$( echo $xyz | cut -f 3 -d ' ' )

    mkdir b$i
    cp b.box b$i/b$i.box
    cp b1e.rea b$i

    cd b$i
    sed -i "5s/.*/-$nex -$ney -$nez/" b$i.box
    genbb b$i > log
    cd ..

  done

  return 0
}

function configure_tests()
{
  export BP_ROOT="$root_dir"/tests/nek5000_bps
  export BENCH_ROOT="$root_dir"

  min_elem=10
  max_elem=12
  min_order=2
  max_order=2

  return 0
}

function build_tests()
{
  # Generate the boxes
  cd $BP_ROOT/boxes
  generate_boxes
  cd "$test_exe_dir"

  # Setup the sin version of the bp
  mkdir sin
  cd sin

  for i in `seq $min_order 1 $max_order`
  do

    mkdir lx$i
    cp -r $BP_ROOT/boxes/b?? $BP_ROOT/SIZE $BP_ROOT/bp1/zsin.usr lx$i/

    # Set lx1 in SIZE file
    sed -i "s/lx1=[0-9]*/lx1=${i}/" lx$i/SIZE

    # Make the executable and copy it into all the
    # box directories
    cd lx$i
    $BP_ROOT/makenek zsin > buildlog
    for j in `seq $min_elem 1 $max_elem`
    do
      cd b$j
      cp ../nek5000 .
      cd ..
    done

    cd ..
  done
}

function run_tests()
{
}

function build_and_run_tests()
{
  configure_tests
  build_tests

  return 0
}

test_required_packages="nek5000"
