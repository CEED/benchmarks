#!/bin/bash

# This file is part of CEED. For more details, see exascaleproject.org.

this_file="${BASH_SOURCE[0]}"
if [[ "${#BASH_ARGV[@]}" -ne "$#" ]]; then
   script_is_sourced="yes"
   exit_cmd=return
else
   script_is_sourced=""
   exit_cmd=exit
fi
test_file=""
config=""
compiler_list=""
build=""
build_list=""
update_packages=""
remove_list=""
run=""
post_process=""
num_proc_build=${num_proc_build:-""}
num_proc_run=${num_proc_run:-""}
num_proc_node=${num_proc_node:-""}
dry_run="" # empty string = NO
start_shell=""
verbose=""
cur_dir="$PWD"

function abspath()
{
   local outvar="$1" path="$2" cur_dir="$PWD"
   cd "$path" && path="$PWD" && cd "$cur_dir" && eval "$outvar=\"$path\""
}

abspath root_dir "$(dirname "$this_file")" || $exit_cmd 1
build_root="$root_dir/builds"
configs_dir="$root_dir/machine-configs"
pkg_sources_dir="$root_dir/package-sources"

if [[ -t 1 ]]; then
   # ANSI color codes
   none=$'\E[0m'
   red=$'\E[0;31m'
   green=$'\E[0;32m'
   yellow=$'\E[0;33m'
   blue=$'\E[0;34m'
   bblue=$'\E[1;34m'
   magenta=$'\E[0;35m'
   cyan=$'\E[0;36m'
   clear="$(tput sgr0)"
fi

help_msg="
$this_file [options]

Options:
   -h|--help             print this usage information and exit
   -c|--config <name>    choose a configuration file
   -m|--compiler \"list\"  choose a compiler, or a list of compilers to use
   -b|--build \"list\"     download and build the listed packages
   -u|--update           update the package sources before (re)building
  -rm|--remove \"list\"    remove packages from the build dir (before building)
   -r|--run <name>       run the tests in the script <name>
   -n|--num-proc \"list\"  total number of MPI tasks to use in the tests
   -p|--proc-node \"list\" number of MPI tasks per node to use in the tests
  -pp|--post-process <name> post process the results using script <name>
   -d|--dry-run          show (but do not run) the commands for the tests
   -s|--shell            execute bash shell commands before running the test
   -v|--verbose          print additional messages
   -x                    enable script tracing with 'set -x'

This script builds and/or runs a set of tests using specified configuration
and compiler.

Example usage:
  $this_file --config vulcan --compiler gcc --build \"metis hypre mfem\"
  $this_file --config vulcan --compiler gcc --run tests/mfem_bps/bp1_v1.sh
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
   [[ -d "$build_root" ]] || mkdir -p "$build_root" || return 1
   output_dir="${build_root}/${short_config}_${compiler}"
   [[ -d "$output_dir" ]] || mkdir -p "$output_dir" || return 1
   OUT_DIR="$output_dir"
   echo "Using OUT_DIR = $OUT_DIR"
   echo
}


function search_file_list()
{
   local out_var="$1" this_file=""
   shift
   for this_file; do
      [[ -e "$this_file" ]] && {
         eval "${out_var}='$this_file'"
         return 0
      }
   done
   eval "${out_var}="
   return 1
}


function add_to_path()
{
   # Usage: add_to_path POS VAR_NAME ["before:"|"after:"|DIR]...
   # POS is "before" or "after". DIR cannot contain ":".
   local out_var="$2" pos="$1" var= item=
   shift 2
   eval var="\${${out_var}}"
   for item; do
      [[ -z "$item" ]] && continue
      case "$item" in
         *:)
            pos="${item%:}"
            continue
            ;;
      esac
      case ":${var}:" in
         *:"$item":*)
            ;;
         ::)
            var="$item"
            ;;
         *)
            if [[ "$pos" = "after" ]]; then
               var="${var}:$item"
            else
               var="${item}:$var"
            fi
            ;;
      esac
   done
   eval "${out_var}=\${var}"
}


function array_union()
{
   # Append an entry to an array, if it is not already in the array.
   # Return false if the entry is already in the array.
   local array_var="$1" entry="$2" idx=""
   eval "local A=(\"\${${array_var}[@]}\")"
   for (( idx=0; idx<"${#A[@]}"; idx++ )); do
      [[ "${A[idx]}" = "$entry" ]] && return 1
   done
   A[idx]="$entry"
   eval "${array_var}=(\"\${A[@]}\")"
   return 0
}


function remove_package()
{
   # Used variables: 'pkg', 'pkg_bld_dir'
   echo "Removing $pkg from OUT_DIR ..."
   rm -rf "${pkg_bld_dir}"{,_build.log,_build_successful}
}


function print_variables()
{
   # Write the variable definitions so that they can be sourced from a file.
   # Usage: print_variables PREFIX [VAR_NAME]...
   local pfx="$1" var_name= var_val= var_list=()
   shift
   for var_name; do
      var_list=("${var_list[@]}" "${pfx}${var_name}")
      eval var_val=\$${var_name}
      var_val=${var_val//\'/"'\\''"}
      printf "${pfx}${var_name}='%s'\n" "$var_val"
   done
   printf "${pfx}VAR_LIST=(%s)\n" "${var_list[*]}"
}


function update_git_package()
{
   # Used variables: 'pkg', 'pkg_src_dir', 'pkg_git_branch', 'pkg_bld_dir'
   if [[ "$update_packages" = "yes" ]]; then
      echo "Updating $pkg ..."
      cd "$pkg_src_dir" && {
         local remote_ref=($(git ls-remote origin $pkg_git_branch))
         local local_ref="$(git rev-parse $pkg_git_branch)"
         if [[ "${remote_ref[0]}" = "$local_ref" ]]; then
            echo "Package $pkg is up to date."
            return 0
         fi
      } && \
      git checkout . && \
      git clean -df && \
      git checkout "$pkg_git_branch" && {
         git pull --ff-only || \
         git checkout -B "$pkg_git_branch" "origin/$pkg_git_branch"
      } || {
         echo "Error updating $pkg. Stop."
         return 1
      }
      cd ..
      : > "${pkg_src_dir}_updated"
   fi
}


function get_package_git_version()
{
   # Used variables: 'pkg_src_dir'
   # Defines the variable: 'pkg_version'
   pkg_version="$(cd "$pkg_sources_dir/$pkg_src_dir" &&
      git describe --long --abbrev=10 --tags 2> /dev/null ||
      git rev-parse HEAD 2> /dev/null)"
   pkg_version="${pkg_version:-(not a git repo?)}"
}


function package_build_is_good()
{
   # Used variables: 'pkg_bld_dir', 'pkg_src_dir', 'pkg', 'pkg_var_prefix'
   local var_name= var_val= cur_var_val= var_list=()
   if [[ -d "${pkg_bld_dir}" && -e "${pkg_bld_dir}_build_successful" ]]; then
      if [[ -e "$pkg_sources_dir/${pkg_src_dir}_updated" ]] && \
         [[ "$pkg_sources_dir/${pkg_src_dir}_updated" -nt \
            "${pkg_bld_dir}_build_successful" ]]; then
         #
         echo "Package $pkg needs to be updated (newer source) ..."
         remove_package
      elif [[ -n "$pkg_var_prefix" ]]; then
         source "${pkg_bld_dir}_build_successful" || {
            echo "Error reading ${pkg_bld_dir}_build_successful ..."
            echo "Package $pkg needs to be updated (read error) ..."
            remove_package
            return 1
         }
         eval "var_list=(\"\${${pkg_var_prefix}VAR_LIST[@]}\")"
         for var_name in "${var_list[@]}"; do
            eval "var_val=\"\${$var_name}\""
            var_name="${var_name#${pkg_var_prefix}}"
            eval "cur_var_val=\"\${$var_name}\""
            [[ -n "$cur_var_val" ]] || cur_var_val="$var_val"
            if [[ "$cur_var_val" != "$var_val" ]]; then
               echo "Package $pkg needs to be updated (modified $var_name) ..."
               remove_package
               return 1
            elif [[ -n "$var_val" && -d "$var_val" && \
                    -e "${var_val}_build_successful" && \
                    "${var_val}_build_successful" -nt \
                    "${pkg_bld_dir}_build_successful" ]]; then
               #
               printf "Package $pkg needs to be updated "
               echo "(newer $(basename "$var_val")) ..."
               remove_package
               return 1
            fi
         done
         return 0
      else
         return 0
      fi
   fi
   return 1
}


function remove_packages()
{
   local _pkg=""
   for _pkg; do
      cd "$root_dir/package-builders"
      if [[ -e "$_pkg.sh" ]]; then
         source "$_pkg.sh" && remove_package || {
            echo "Error removing package \"$_pkg\". Stop."
            return 1
         }
      else
         echo "Package \"$_pkg\" does not exist. Stop."
         return 1
      fi
   done
}


function build_packages()
{
   local _pkg=
   mkdir -p "$pkg_sources_dir"
   for _pkg; do
      cd "$root_dir/package-builders"
      if [[ -e "$_pkg.sh" ]]; then
         unset -v pkg_version
         unset -v pkg_var_prefix
         unset -f build_package
         source "$_pkg.sh" && build_package || {
            echo "Error building package \"$_pkg\". Stop."
            return 1
         }
         echo "$_pkg version: $pkg_version"
      else
         echo "Package \"$_pkg\" does not exist. Stop."
         return 1
      fi
   done
   cd "$cur_dir"
}


function compose_mpi_run_command()
{
   mpi_run="${MPIEXEC:-mpirun} ${MPIEXEC_OPTS}"
   mpi_run+=" ${MPIEXEC_NP:--np} ${num_proc_run} $bind_sh"
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
         return 1
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
   [[ "$verbose" = "yes" ]] || return 0
   echo "Compilers:"
   echo "----------"
   which $MPICC
   echo
   $MPICC $mpi_info_flag
   echo
   which $MPICXX
   echo
   $MPICXX $mpi_info_flag
   echo
}


### Process command line parameters

while [ $# -gt 0 ]; do

case "$1" in
   -h|--help)
      # Echo usage information
      echo "$help_msg"
      print_configs
      $exit_cmd
      ;;
   -c|--config)
      shift
      [ $# -gt 0 ] || { echo "Missing <name> in --config <name>"; $exit_cmd 1; }
      config="$1"
      [ -r "$config" ] || {
         config="$configs_dir/${config}.sh"
         [ -r "$config" ] || {
            echo "Configuration file not found: '$1' / '$config'"; $exit_cmd 1
         }
      }
      ;;
   -m|--compiler)
      shift
      [ $# -gt 0 ] || {
      echo "Missing \"list\" in --compiler \"list\""; $exit_cmd 1; }
      compiler_list="$1"
      ;;
   -b|--build)
      build=on
      shift
      [ $# -gt 0 ] || {
      echo "Missing \"list\" in --build \"list\""; $exit_cmd 1; }
      build_list="$1"
      ;;
   -u|--update)
      update_packages="yes"
      ;;
   -rm|--remove)
      shift
      [ $# -gt 0 ] || {
      echo "Missing \"list\" in --remove \"list\""; $exit_cmd 1; }
      remove_list="$1"
      ;;
   -r|--run)
      run=on
      shift
      [ $# -gt 0 ] || { echo "Missing <name> in --run <name>"; $exit_cmd 1; }
      test_file="$1"
      [[ -r "$test_file" ]] || {
         echo "Test script not found: '$1'"; $exit_cmd 1
      }
      ;;
   -n|--num-proc)
      shift
      [ $# -gt 0 ] || {
      echo "Missing \"list\" in --num-proc \"list\""; $exit_cmd 1; }
      num_proc_run="$1"
      ;;
   -p|--proc-node)
      shift
      [ $# -gt 0 ] || {
      echo "Missing \"list\" in --proc-node \"list\""; $exit_cmd 1; }
      num_proc_node="$1"
      ;;
   -pp|--post-process)
      post_process=on
      shift
      [ $# -gt 0 ] || { echo "Missing <name> in --post-process <name>"; $exit_cmd 1; }
      pp_file="$1"
      [[ -r "$pp_file" ]] || {
         echo "Post process script not found: '$1'"; $exit_cmd 1
      }
      ;;
   -d|--dry-run)
      dry_run="quoted_echo"
      ;;
   -s|--shell)
      start_shell="yes"
      ;;
   -v|--verbose)
      verbose="yes"
      ;;
   -x)
      set -x
      ;;
   *)
      echo "Unknown option: '$1'"
      $exit_cmd 1
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
   $exit_cmd 1
}

echo "Reading configuration $config ..."
. "$config" || $exit_cmd 1

[ -z "$valid_compilers" ] && {
   echo "Invalid configuration file: $config"
   $exit_cmd 1
}

abspath config_dir "$(dirname "$config")" || $exit_cmd 1
short_config="$(basename "$config")"
config="${config_dir}/${short_config}"
short_config="${short_config#config_}"
short_config="${short_config%.sh}"

[ -z "$compiler_list" ] && {
   echo
   echo "Choose compiler(s) with -m \"list\""
   echo
   echo "Available compilers: $valid_compilers"
   echo
   $exit_cmd 1
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
   $exit_cmd 1
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
   $exit_cmd 1
}

echo "Setting up compiler $compiler ..."
mpi_info_flag=-v
setup_"$compiler"
echo
show_compilers


set_build_dirs || $exit_cmd 1


### Remove packages

[[ -n "$remove_list" ]] && {
   remove_packages $remove_list || $exit_cmd 1
   echo
}


### Build packages

[[ -n "$build" ]] && {
   num_proc_build=${num_proc_build:-4}
   echo "Building packages using $num_proc_build processors."

   build_packages $build_list || $exit_cmd 1
   echo
}


### Run the tests (building and running $test_file)

[ -n "$run" ] && {

cd "$cur_dir"
abspath test_dir "$(dirname "$test_file")" || $exit_cmd 1
test_basename="$(basename "$test_file")"
test_file="${test_dir}/${test_basename}"

[[ "$verbose" = "yes" ]] && {
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
}

test_up_dir="${test_dir#$root_dir/tests}"
test_exe_dir="$OUT_DIR/${test_up_dir#/}"
echo "Creating test executables directory: OUT_DIR/${test_up_dir#/}"
$dry_run mkdir -p "$test_exe_dir" || $exit_cmd 1

trap 'printf "\nScript interrupted.\n"; '$exit_cmd' 33' INT

## Source the test script file.
echo "Reading test file: $test_file"
echo
test_required_packages=""
. "$test_file" || $exit_cmd 1

## Build any packages required by the test
echo "Packages required by the test: $test_required_packages"
build_packages $test_required_packages || $exit_cmd 1
echo

## Loop over the number-of-processors list.
for (( num_proc_idx = 0; num_proc_idx < num_proc_list_size; num_proc_idx++ ))
do

num_proc_run="${num_proc_list[$num_proc_idx]}"
num_proc_node="${num_proc_node_list[$num_proc_idx]}"

set_num_nodes || $exit_cmd 1

if [[ "$start_shell" = "yes" ]]; then
   if [[ ! -t 1 ]]; then
      echo "Standard output is not a terminal. Stop."
      $exit_cmd 1
   fi
   echo "Reading shell commands, type 'c' to continue, 'exit' to stop ..."
   echo
   cd "$cur_dir"
   set -o emacs
   PS1='$ '
   [[ -r $HOME/.bashrc ]] && source $HOME/.bashrc
   HISTFILE="$root_dir/.bash_history"
   history -c
   history -r
   # bind '"\\C-i": menu-complete'
   alias c='break'
   while cwd="$PWD/" cwd="${cwd#${root_dir}/}" cwd="${cwd%/}" \
         prompt="[${cyan}benchmarks$none:$blue$cwd$clear]\$ " && \
         read -p "$prompt" -e line; do
      history -s "$line"
      history -w
      shopt -q -s expand_aliases
      eval "$line"
      shopt -q -u expand_aliases
   done
   [[ "${#line}" -eq 0 ]] && { echo; $exit_cmd 0; }
   shopt -q -u expand_aliases
   echo "Continuing ..."
fi

build_and_run_tests
echo

done ## End of loop over processor numbers

trap - INT

} ## run is on

### Post process the results

[[ -n "$post_process" ]] && {

. "$pp_file" || $exit_cmd 1

abspath test_dir "$(dirname "$pp_file")" || $exit_cmd 1
test_up_dir="${test_dir#$root_dir/tests}"
test_exe_dir="$OUT_DIR/${test_up_dir#/}"

postprocess

} ## post-process on 

) || $exit_cmd 1
done ## Loop over $compiler_list


$exit_cmd 0
