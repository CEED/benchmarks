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
}

function generate_boxes()
{
  cd $BP_ROOT/boxes
  # Run thorugh the box sizes
  for i in `seq $min_elem 1 $max_elem`
  do
    # Generate the boxes only if they have not
    # been generated before.
    if [[ ! -f b$i/b$i.rea ]]; then
      # Set the number of elements in box file.
      xyz=$(xyz $i)
      nex=$( echo $xyz | cut -f 1 -d ' ' )
      ney=$( echo $xyz | cut -f 2 -d ' ' )
      nez=$( echo $xyz | cut -f 3 -d ' ' )

      mkdir -p b$i
      sed "5s/.*/-$nex -$ney -$nez/" b.box > b$i/b$i.box
      cp b1e.rea b$i

      cd b$i
      genbb b$i &> log
      cd ..
    fi
  done
}


function configure_tests()
{
  export BP_ROOT="$root_dir"/tests/nek5000_bps

  min_elem=6
  max_elem=122
  min_order=3
  max_order=4
}

function build_tests()
{
  cd "$test_exe_dir"

  # Setup the sin version of the bp
  mkdir -p $1 
  cd $1 

  # Export variables needed by the 'makenek' script.
  local CFLAGS_orig="$CFLAGS" FFLAGS_orig="$FFLAGS"
  CFLAGS="${CFLAGS//:/\\:}"
  FFLAGS="${FFLAGS//:/\\:}"
  PPLIST="$NEK5K_EXTRA_PPLIST"
  # export NEK5K_DIR CFLAGS FFLAGS MPIF77 MPICC PPLIST
  export NEK5K_DIR MPIF77 MPICC PPLIST CFLAGS FFLAGS

  for i in `seq $min_order 1 $max_order`
  do
    # Only build nek5000 if it is not built
    # already.
    if [[ ! -e lx$i ]]; then
      mkdir -p lx$i
      cp -r $BP_ROOT/boxes/b?* $BP_ROOT/bp1/$1.usr lx$i/

      # Set lx1 in SIZE file
      sed "s/lx1=[0-9]*/lx1=${i}/" $BP_ROOT/SIZE > lx$i/SIZE

      # Make the executable and copy it into all the
      # box directories
      cd lx$i
      echo "Building the $1 tests in directory $PWD ..."
      $BP_ROOT/makenek $1 &> buildlog
      if [[ ! -e nek5000 ]]; then
        echo "Error building the test, see 'buildlog' for details. Stop."
        CFLAGS="${CFLAGS_orig}"
        FFLAGS="${FFLAGS_orig}"
        return 1
      fi
      for j in `seq $min_elem 1 $max_elem`
      do
        cp ./nek5000 b$j/
      done

      cd ..
    fi
  done

  cd ..

  cp -r $1 $2
  cd $2

  for i in `seq $min_order 1 $max_order`
  do
    cp -r $BP_ROOT/bp1/$2.usr lx$i/
    cd lx$i
    rm nek5000 $1.usr

    echo "Building the $2 tests in directory $PWD ..."
    $BP_ROOT/makenek $2 &> buildlog
    if [[ ! -e nek5000 ]]; then
      echo "Error building the test, see 'buildlog' for details. Stop."
      CFLAGS="${CFLAGS_orig}"
      FFLAGS="${FFLAGS_orig}"
      return 1
    fi

    for j in `seq $min_elem 1 $max_elem`
    do
      cp ./nek5000 b$j/
    done

    cd ..
  done

  CFLAGS="${CFLAGS_orig}"
  FFLAGS="${FFLAGS_orig}"
}

function nekmpi()
{
  cp $BP_ROOT/"submit.sh" .

  # This next check is important only for machines without virtual memory like
  # Blue Gene/Q.
  local size_out=($(size ./nek5000))
  local proc_size=${size_out[9]}
  if [[ -n "$node_virt_mem_lim" ]] && \
     (( proc_size * num_proc_node > node_virt_mem_lim * 1024*1024*1024 )); then
     echo "The total section size in ./nek5000 is too big: $proc_size ..."
     echo " ... in directory: $PWD"
     return 1
  fi

  ./submit.sh $1 $2
}

function run_tests()
{
  cd "$test_exe_dir"

  set_mpi_options
  local mpi_run="${MPIEXEC:-mpirun} $MPIEXEC_OPTS"
  export mpi_run="$mpi_run ${MPIEXEC_NP:--np} $num_proc_run $bind_sh"

  cd $1 

  for i in `seq $min_order 1 $max_order`
  do
    cd lx$i
    for j in `seq $min_elem 1 $max_elem`
    do
      cd b$j
      nekmpi b$j $num_proc_run
      cd ..
    done

    cd ..
  done

  cd ..
}

function build_and_run_tests()
{
  echo 'Setting up the tests ...'
  $dry_run configure_tests
  echo "Generating the box meshes ..."
  $dry_run generate_boxes
  echo 'Buiding the sin and w tests ...'
  $dry_run build_tests zsin zw || return 1
  echo 'Running the sin tests ...'
  $dry_run run_tests zsin
  echo 'Running the w tests ...'
  $dry_run run_tests zw
}

test_required_packages="nek5000"
