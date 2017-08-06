#!/bin/bash

echo $2

if [ -z $3 ]
then
	echo 'Please enter the command as: check.sh [muse/vimos] [galaxy] [opt_MCdir] [opptional dir]'
	exit	
fi

cd /Data/$1/analysis/$2/$3/MC/$4
max=$( awk '{print $3}' /Data/$1/analysis/$2/$3/setup/voronoi_2d_binning_output.txt | sort -n | tail -1 )




i=0
while [[ $i -lt $max+1 ]]
do 
    ls ${i}.dat &> /dev/null || echo "     "$i / $max
    let i="i+1"
done
# done
