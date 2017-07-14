# Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
# the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
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


source "$root_dir"/tests/nek5000_bps/make-boxes.sh


function configure_tests()
{
  export BP_ROOT="$root_dir"/tests/nek5000_bps

  # Settings for 1 node:

  # {min,max}_elem are the exponents for the number of elements, base is 2
  min_elem=1
  max_elem=21
  # the "order" here is actually number of 1D points, i.e. p+1, not p
  min_order=2
  max_order=9
  # the number of points is computed as num_elements*(p+1)**3
  max_points=3000000

  while (( 2**min_elem < num_proc_node )); do
     ((min_elem=min_elem+1))
  done

  set_max_elem_order "$min_order"
  max_elem="$max_elem_order"

  # Settings for more than 1 node:

  local n=$num_nodes
  while (( n >= 2 )); do
     ((min_elem=min_elem+1))
     ((max_elem=max_elem+1))
     ((max_points=2*max_points))
     ((n=n/2))
  done

  # Make sure that we do not exceed 2^21 limit
  if [[ "$max_elem" -gt 21 ]]; then
    max_elem=21
  fi
}

function set_max_elem_order()
{
  max_elem_order="$min_elem"
  local pp1="$1" s=
  for ((s = min_elem; s <= max_elem; s++)); do
    local npts=$(( 2**s * (pp1-1)**3 ))
    (( npts > max_points )) && break
    max_elem_order="$s"
  done
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
    newbuild=false
    if [[ ! -e "lx$i" ]]; then
      mkdir -p lx$i
      newbuild=true
    fi

    set_max_elem_order "$i"
    
    local lelg=$(( 2**max_elem_order )) 
    local lelt=$(( lelg/num_proc_run ))
    sed -e "s/lelt=[0-9]*/lelt=${lelt}/" \
        -e "s/lp=[0-9]*/lp=${num_proc_run}/" \
        -e "s/lelg=[0-9]*/lelg=${lelg}/" \
        $BP_ROOT/SIZE > SIZE

    cp $BP_ROOT/bp1/$1.usr lx$i/

    # Set lx1 in SIZE file
    sed "s/lx1=[0-9]*/lx1=${i}/" SIZE > lx$i/SIZE.new
    if [[ "$newbuild" == "true" ]]; then
      mv "lx$i"/SIZE.new "lx$i"/SIZE
    elif [[ ! -f "lx$i"/SIZE ]]; then
      mv "lx$i"/SIZE.new "lx$i"/SIZE
      newbuild=true
    else
      if ! cmp -s "lx$i"/SIZE "lx$i"/SIZE.new; then
        newbuild=true
        rm "lx$i"/SIZE
        mv "lx$i"/SIZE.new "lx$i"/SIZE
      fi
    fi

    # Make the executable and copy it into all the
    # box directories
    cd lx$i
    echo "Building the $1 tests in directory $PWD ..."
    if [[ "$newbuild" = false ]]; then
      echo "Reusing the existing build ..."
    fi
    $BP_ROOT/makenek $1 &> buildlog
    if [[ ! -e nek5000 ]]; then
      echo "Error building the test, see 'buildlog' for details. Stop."
      CFLAGS="${CFLAGS_orig}"
      FFLAGS="${FFLAGS_orig}"
      return 1
    fi

    for j in `seq $min_elem 1 $max_elem_order`
    do
      cp -r $BP_ROOT/boxes/b$j .
      cp ./nek5000 b$j/
    done

    cd .. ## lx$i
  done

  cd ..

  if [[ "$2" != "" ]]; then
    cp -r $1 $2
    cd $2

    for i in `seq $min_order 1 $max_order`
    do
      cp -r $BP_ROOT/bp1/$2.usr lx$i/
      cd lx$i
      rm nek5000 $1.usr > /dev/null 2>&1

      echo "Building the $2 tests in directory $PWD ..."
      $BP_ROOT/makenek $2 &> buildlog
      if [[ ! -e nek5000 ]]; then
        echo "Error building the test, see 'buildlog' for details. Stop."
        CFLAGS="${CFLAGS_orig}"
        FFLAGS="${FFLAGS_orig}"
        return 1
      fi

      set_max_elem_order "$i"

      for j in `seq $min_elem 1 $max_elem_order`
      do
        cp ./nek5000 b$j/
      done

      cd ..
    done
  fi

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
  [[ -z "$dry_run" ]] && cd "$test_exe_dir"

  set_mpi_options
  export mpi_run

  [[ -z "$dry_run" ]] && cd $1

  for i in `seq $min_order 1 $max_order`
  do
    [[ -z "$dry_run" ]] && cd lx$i
    set_max_elem_order "$i"
    for j in `seq $min_elem 1 $max_elem_order`
    do
      local npts=$(( ((i-1)**3) * (2**j) ))
      echo
      printf "Running order $((i-1)), with number of elements 2^$j;"
      echo " number of points is $npts."

      [[ -z "$dry_run" ]] && cd b$j
      $dry_run nekmpi b$j $num_proc_run
      [[ -z "$dry_run" ]] && cd ..
    done

    [[ -z "$dry_run" ]] && cd ..
  done

  [[ -z "$dry_run" ]] && cd ..
}

function build_and_run_tests()
{
  echo 'Setting up the tests ...'
  configure_tests
  echo "Generating the box meshes ..."
  $dry_run generate_boxes || return 1
  echo 'Buiding the sin and w tests ...'
  $dry_run build_tests zsin || return 1
#  $dry_run build_tests zsin zw || return 1
  echo 'Running the sin tests ...'
  run_tests zsin
  # W tests are commented as there is no diskquota
  # in vulcan to run both the tests
#  echo 'Running the w tests ...'
#  $dry_run run_tests zw
}

test_required_packages="nek5000"
