#!/bin/bash

echo $1

if [ -z "$3" ]
then
    cd /Data/vimos/analysis/$1/gas_MC/
    max=$( awk '{print $3}' ../voronoi_2d_binning_output_$2.txt | sort -n | tail -1 )
else
    cd /Data/vimos/analysis/$1/gas_MC/$3/
    max=$( awk '{print $3}' ../../voronoi_2d_binning_output_$2.txt | sort -n | tail -1 )
fi

i=0
while [[ $i -lt $max ]]
do 
    ls ${i}.dat &> /dev/null || echo "     "$i / $max
    let i="i+1"
done
# done
