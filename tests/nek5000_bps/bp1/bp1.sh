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

  min_elem=10
  max_elem=12
  min_order=2
  max_order=2
}

function build_tests()
{
  # Generate the boxes
  cd $BP_ROOT/boxes
  echo "Generating the box meshes ..."
  generate_boxes
  cd "$test_exe_dir"

  # Setup the sin version of the bp
  mkdir -p sin
  cd sin

  # Export variables needed by the 'makenek' script.
  # CFLAGS="${CFLAGS//:/\\:}"
  # FFLAGS="${FFLAGS//:/\\:}"
  PPLIST="$NEK5K_EXTRA_PPLIST"
  # export NEK5K_DIR CFLAGS FFLAGS MPIF77 MPICC PPLIST
  export NEK5K_DIR MPIF77 MPICC PPLIST

  for i in `seq $min_order 1 $max_order`
  do
    # Only build nek5000 if it is not built
    # already.
    if [[ ! -e lx$i ]]; then 
      mkdir -p lx$i
      cp -r $BP_ROOT/boxes/b?? $BP_ROOT/bp1/zsin.usr lx$i/

      # Set lx1 in SIZE file
      sed "s/lx1=[0-9]*/lx1=${i}/" $BP_ROOT/SIZE > lx$i/SIZE

      # Make the executable and copy it into all the
      # box directories
      cd lx$i
      echo "Building the test in directory $PWD ..."
      $BP_ROOT/makenek zsin &> buildlog
      if [[ ! -e nek5000 ]]; then
        echo "Error building the test, see 'buildlog' for details. Stop."
        return 1
      fi
      for j in `seq $min_elem 1 $max_elem`
      do
        cd b$j
        cp ../nek5000 .
        cd ..
      done

      cd ..
    fi
  done

  cd ..
}

function nekmpi()
{
  cp $BP_ROOT/"submit.sh" .

  if [[ "$short_config" = "vulcan" ]]; then
    sbatch ./submit.sh $1 $2
  elif [[ "$short_config" = "linux" || "$short_config" = "mac" ]]; then
    ./submit.sh $1 $2 
  fi
}

function run_tests()
{
  cd "$test_exe_dir"

  set_mpi_options
  local mpi_run="${MPIEXEC:-mpirun} $MPIEXEC_OPTS"
  export mpi_run="$mpi_run ${MPIEXEC_NP:--np} $num_proc_run $bind_sh"

  cd sin

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
  echo 'Buiding the tests ...'
  $dry_run build_tests || return 1
  echo 'Running the tests ...'
  $dry_run run_tests
}

test_required_packages="nek5000"
