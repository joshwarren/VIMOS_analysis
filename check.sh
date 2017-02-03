#!/bin/bash

echo $1

if [ $2 = kin ]
then 
	cd /Data/vimos/analysis/$1/gas_MC/$3/
elif [ $2 = pop ]
then
	cd /Data/vimos/analysis/$1/pop_MC/$3/
fi

max=$( awk '{print $3}' /Data/vimos/analysis/$1/voronoi_2d_binning_output_$2.txt | sort -n | tail -1 )


i=0
while [[ $i -lt $max ]]
do 
    ls ${i}.dat &> /dev/null || echo "     "$i / $max
    let i="i+1"
done
# done
