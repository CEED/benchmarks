#!/bin/bash

# This file is part of CEED. For more details, see exascaleproject.org.

test_file=""
config=""
compiler_list=""
build=""
run=""
num_proc_build=""
num_proc_run=""
num_proc_node=""
dry_run="" # empty string = NO
cur_dir="$PWD"
root_dir="$(dirname "$0")"
cd "$root_dir" && root_dir="$PWD" && cd "$cur_dir" || exit 1
build_root="$root_dir/builds"
configs_dir="$root_dir/machine-configs"
pkg_sources_dir="$root_dir/package-sources"

help_msg="
$0 [options]

Options:
   -h|--help             print this usage information and exit
   -c|--config <name>    choose a configuration file
   -m|--compiler \"list\"  choose a compiler, or a list of compilers to use
   -b|--build \"list\"     download and build the listed packages
   -r|--run <name>       run the tests in the script <name>
   -n|--num-proc \"list\"  total number of MPI tasks to use in the tests
   -p|--proc-node \"list\" number of MPI tasks per node to use in the tests
   -d|--dry-run          show (but do not run) the commands for the tests

This script builds and/or runs a set of tests using specified configuration
and compiler.

Example usage:
  $0 --config vulcan --compiler gcc --build \"metis hypre mfem\"
  $0 --config vulcan --compiler gcc --run tests/mfem_bps/bp1_v1.sh
"

function print_configs()
{
   configs=""
   for conf_file in "$configs_dir"/*.sh; do
      [ -r "$conf_file" ] && {
         conf="${conf_file#$configs_dir/}"
         conf="${conf%.sh}"
         [ -n "$conf" ] && {
            [ -z "$configs" ] && configs="$conf" || \
            configs="$configs $conf"
         }
      }
   done
   [ -z "$configs" ] && configs="(no configuration files found)"
   echo "Available configurations: $configs"
   echo
}


function set_build_dirs()
{
   # Setup separate build directories inside $build_root based on $config
   # and $compiler.
   [[ -d "$build_root" ]] || mkdir -p "$build_root" || exit 1
   output_dir="${build_root}/${short_config}_${compiler}"
   [[ -d "$output_dir" ]] || mkdir -p "$output_dir" || exit 1
   cd "$output_dir" && OUT_DIR="$PWD" && cd "$cur_dir" || exit 1
   echo "Using OUT_DIR = $OUT_DIR"
   echo
}


function build_packages()
{
   mkdir -p "$pkg_sources_dir"
   for pkg; do
      cd "$root_dir/package-builders"
      if [[ -e "$pkg.sh" ]]; then
         . "$pkg.sh" || {
            echo "Error building package \"$pkg\". Stop."
            exit 1
         }
      else
         echo "Package \"$pkg\" does not exist. Stop."
         exit 1
      fi
   done
   cd "$cur_dir"
}


function check_memory_req()
{
   local total_mem=""
   if [[ -n "$num_nodes" && -n "$memory_per_node" && \
         -n "$total_memory_required" ]]; then
      ((total_mem = memory_per_node * num_nodes))
      # echo "Total memory available: $total_mem GiB"
      # echo "Total memory required : $total_memory_required GiB"
      if [[ "$total_memory_required" -gt "$total_mem" ]]; then
         printf " *** Insufficient total memory: $total_mem GiB, "
         printf "this test requires: $total_memory_required GiB. "
         echo "Skipping test."
         return 1
      fi
   else
      echo " *** Warning: unable to check memory requirement."
   fi
   return 0
}


function quoted_echo()
{
   local arg= string=
   for arg; do
      if [[ -z "${arg##* *}" ]]; then
         string+=" \"${arg//\"/\\\"}\""
      else
         string+=" $arg"
      fi
   done
   printf "%s\n" "${string# }"
}


function set_num_nodes()
{
   if [[ -n "$num_proc_node" ]]; then
      ((num_proc_run % num_proc_node != 0)) && {
         echo "The total number of tasks ($num_proc_run) must be a multiple of"
         echo "the number of tasks per node ($num_proc_node). Stop."
         exit 1
      }
      ((num_nodes = num_proc_run / num_proc_node))
   else
      num_proc_node="unknown number of"
      num_nodes=""
   fi
   echo "Running the tests using a total of $num_proc_run MPI tasks ..."
   echo "... with $num_proc_node tasks per node ..."
   echo
}


function show_compilers()
{
   echo "Compilers:"
   echo "----------"
   which $mpi_cc
   echo
   $mpi_cc $mpi_info_flag
   echo
   which $mpi_cxx
   echo
   $mpi_cxx $mpi_info_flag
   echo
}


### Process command line parameters

while [ $# -gt 0 ]; do

case "$1" in
   -h|--help)
      # Echo usage information
      echo "$help_msg"
      print_configs
      exit
      ;;
   -c|--config)
      shift
      [ $# -gt 0 ] || { echo "Missing <name> in --config <name>"; exit 1; }
      config="$1"
      [ -r "$config" ] || {
         config="$configs_dir/${config}.sh"
         [ -r "$config" ] || {
            echo "Configuration file not found: '$1' / '$config'"; exit 1
         }
      }
      ;;
   -m|--compiler)
      shift
      [ $# -gt 0 ] || {
      echo "Missing \"list\" in --compiler \"list\""; exit 1; }
      compiler_list="$1"
      ;;
   -b|--build)
      build=on
      shift
      [ $# -gt 0 ] || {
      echo "Missing \"list\" in --build \"list\""; exit 1; }
      build_list="$1"
      ;;
   -r|--run)
      run=on
      shift
      [ $# -gt 0 ] || { echo "Missing <name> in --run <name>"; exit 1; }
      test_file="$1"
      [[ -r "$test_file" ]] || {
         echo "Test script not found: '$1'"; exit 1
      }
      ;;
   -n|--num-proc)
      shift
      [ $# -gt 0 ] || {
      echo "Missing \"list\" in --num_proc \"list\""; exit 1; }
      num_proc_run="$1"
      ;;
   -p|--proc-node)
      shift
      [ $# -gt 0 ] || {
      echo "Missing \"list\" in --proc-node \"list\""; exit 1; }
      num_proc_node="$1"
      ;;
   -d|--dry-run)
      dry_run="quoted_echo"
      ;;
   *)
      echo "Unknown option: '$1'"
      exit 1
      ;;
esac

shift
done # while ...
# Done processing command line parameters


### Read configuration file

[ -z "$config" ] && {
   echo
   echo "Choose a configuration with -c <name>, or use -h for help"
   echo
   print_configs
   exit 1
}

echo "Reading configuration $config ..."
. "$config" || exit 1

[ -z "$valid_compilers" ] && {
   echo "Invalid configuration file: $config"
   exit 1
}

config_dir="$(dirname "$config")"
cd "${config_dir}" && config_dir="$PWD" && cd "$cur_dir" || exit 1
short_config="$(basename "$config")"
config="${config_dir}/${short_config}"
short_config="${short_config#config_}"
short_config="${short_config%.sh}"

test_dir="$(dirname "$test_file")"
cd "${test_dir}" && test_dir="$PWD" && cd "$cur_dir" || exit 1
test_basename="$(basename "$test_file")"
test_file="${test_dir}/${test_basename}"

[ -z "$compiler_list" ] && {
   echo
   echo "Choose compiler(s) with -m \"list\""
   echo
   echo "Available compilers: $valid_compilers"
   echo
   exit 1
}

num_proc_list=(${num_proc_run:-4})
num_proc_list_size=${#num_proc_list[@]}
num_proc_node_list=(${num_proc_node:-4})
num_proc_node_list_size=${#num_proc_node_list[@]}
(( num_proc_list_size != num_proc_node_list_size )) && {
   echo "
The size of the number-of-processors list (option --num-proc) must be the same
as the size of the number-of-processors-per-node list (option --proc-node)."
   echo
   exit 1
}


### Loop over compilers

for compiler in $compiler_list; do
(  ## Run each compiler in its own environment

### Setup the environment based on $compiler

# Check if $compiler is valid
valid_comp_pat=" ${valid_compilers} "
[ -n "${valid_comp_pat##* $compiler *}" ] && {
   echo
   echo "Invalid compiler: '$compiler'"
   echo
   echo "Available compilers: $valid_compilers"
   echo
   exit 1
}

echo "Setting up compiler $compiler ..."
mpi_info_flag=-v
setup_"$compiler"
echo
show_compilers


### Build packages

set_build_dirs

[[ -n "$build" ]] && {
   num_proc_build=${num_proc_build:-4}
   echo "Building packages using $num_proc_build processors."

   build_packages $build_list
   echo
}


### Run the tests (building and running $test_file)

[ -n "$run" ] && {

echo "Config file, $(basename "$config"):"
echo "------------------------------------------------"
cat $config
echo "------------------------------------------------"
echo

echo "Test problem file, $test_basename:"
echo "------------------------------------------------"
cat $test_file
echo "------------------------------------------------"
echo

test_up_dir="$(basename "$test_dir")"
test_exe_dir="$OUT_DIR/$test_up_dir"
echo "Creating test executables directory: OUT_DIR/$test_up_dir"
$dry_run mkdir -p "$test_exe_dir" || exit 1

trap 'printf "\nScript interrupted.\n"; exit 33' INT

## Source the test script file.
test_required_packages=""
. "$test_file" || exit 1

## Build any packages required by the test
build_packages $test_required_packages
echo

## Loop over the number-of-processors list.
for (( num_proc_idx = 0; num_proc_idx < num_proc_list_size; num_proc_idx++ ))
do

num_proc_run="${num_proc_list[$num_proc_idx]}"
num_proc_node="${num_proc_node_list[$num_proc_idx]}"

set_num_nodes

build_and_run_tests
echo

done ## End of loop over processor numbers

trap - INT

} ## run is on

) || exit 1
done ## Loop over $compiler_list


exit 0
