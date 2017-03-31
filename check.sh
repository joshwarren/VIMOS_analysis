#!/bin/bash

echo $1

if [ -z $2 ]
then
	echo 'Please enter the command as: check.sh [galaxy] [opt] [opptional dir]'
	exit	
fi

if [ $2 = kin ]
then 
	cd /Data/vimos/analysis/$1/gas_MC/$3/
	binsfile=/Data/vimos/analysis/$1/voronoi_2d_binning_output_$2.txt
elif [ $2 = abs ]
then
	cd /Data/vimos/analysis/$1/pop_MC/$3/
	binsfile=/Data/vimos/analysis/$1/voronoi_2d_binning_output_pop.txt

elif [ $2 = pop ]
then
	cd /Data/vimos/analysis/$1/pop/$3/
	binsfile=/Data/vimos/analysis/$1/voronoi_2d_binning_output_$2.txt
fi

max=$( awk '{print $3}' $binsfile | sort -n | tail -1 )


i=0
while [[ $i -lt $max+1 ]]
do 
    ls ${i}.dat &> /dev/null || echo "     "$i / $max
    # ls ${i}.dat &> /dev/null || scp warrenj@glamdring.physics.ox.ac.uk:analysis/$1/pop/$i.dat /Data/vimos/analysis/$1/pop/
    # ls ${i}.dat &> /dev/null || scp warrenj@glamdring.physics.ox.ac.uk:analysis/$1/pop/distribution/$i.dat /Data/vimos/analysis/$1/pop/distribution/
    # ls ${i}.dat &> /dev/null || scp warrenj@glamdring.physics.ox.ac.uk:analysis/$1/pop/plots/$i.png /Data/vimos/analysis/$1/pop/plots/

    let i="i+1"
done
# done
