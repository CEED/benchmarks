#!/bin/bash -x

# Generate boxes for the bakeoff problems

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

# Run thorugh the box sizes
for i in `seq 10 1 12`                                           # nelt
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
sed -i "5s/.*/-$nex -$ney -$nez/" 'b'$i.box
genbb 'b'$i
cd ..

done
