#!/bin/bash

# for gal in  ngc3557 ic1459 ic1531 ic4296 ngc0612 ngc1399 ngc3100 ngc7075 \
#     pks0718-34 eso443-g024
for gal in  eso443-g024
do
    echo $gal

    if [ -n "$2" ]
    then
        cd /Data/vimos/analysis/$gal/gas_MC/
        max=$( awk '{print $3}' ../voronoi_2d_binning_output_$1.txt | sort -n | tail -1 )
    else
        cd /Data/vimos/analysis/$gal/gas_MC/$2/
        max=$( awk '{print $3}' ../../voronoi_2d_binning_output_$1.txt | sort -n | tail -1 )
    fi

    i=0
    while [[ $i -lt $max ]]
    do 
        ls ${i}.dat &> /dev/null || echo "     "$i / $max
        let i="i+1"
    done
done
