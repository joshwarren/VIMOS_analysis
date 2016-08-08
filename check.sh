#!/bin/bash


for gal in  ngc3557 ic1459 ic1531 ic4296 ngc0612 ngc1399 ngc3100 ngc7075 \
    pks0718-34 eso443-g024
do
    echo $gal
    cd /cygdrive/x/Data/vimos/analysis/$gal/gas_MC/
    max=$( awk '{print $3}' ../voronoi_2d_binning_output.txt | sort -n | tail -1 )
    i=0
    while [[ $i -lt $max ]]
    do 
        ls ${i}.dat &> /dev/null || echo "     "$i / $max
        let i="i+1"
    done
done
