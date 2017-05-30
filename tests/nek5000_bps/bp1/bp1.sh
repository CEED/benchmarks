#!/bin/bash

export BP_ROOT=`pwd`/..
export BENCH_ROOT="$BP_ROOT/../.."

# Generate the boxes
cd $BP_ROOT/boxes
./boxes.sh
cd $BP_ROOT/bp1

# Setup the sin version of the bp
mkdir sin
cd sin

for i in `seq 2 1 4`
do

mkdir lx$i
cp -r $BP_ROOT/boxes/b?? $BP_ROOT/SIZE ../zsin.usr lx$i/

# Set lx1 in SIZE file
sed -i "s/lx1=[0-9]*/lx1=${i}/" lx$i/SIZE

# Make the executable
cd lx$i
$BP_ROOT/makenek zsin
cd ..

done
