#/bin/bash

mkdir w
cp zw SIZE makenek submit.sh w/
cp -r b*  w/

mkdir sin
cp zsin SIZE makenek submit.sh sin/
cp -r b*  w/

cd w
for lx1 in `seq 2 1 11`
do

# Set lx1 in SIZE file
sed -i "s/lx1=[0-9]*/lx1=${lx1}" SIZE

# Make the executable
./makenek zw

for n in `seq 14 1 21`
do
cd 'b'${n}
cp ../nek5000 .
cp ../submit.sh .
qsub -n 512 ./submit.sh 'b'${n} 32
sleep 5
cd ..
done

done

cd ..
cd sin
for lx1 in `seq 2 1 11`
do

# Set lx1 in SIZE file
sed -i "s/lx1=[0-9]*/lx1=${lx1}" SIZE

# Make the executable
./makenek zsin

for n in `seq 14 1 21`
do
cd 'b'${n}
cp ../nek5000 .
cp ../submit.sh .
qsub -n 512 ./submit.sh 'b'${n} 32
sleep 5
cd ..
done

done
