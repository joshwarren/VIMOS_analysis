#!/bin/bash



imhead VIMOS.2012-06-18T09\:29\:07.127.fits | grep "HIERARCH ESO DPR TYPE" | cut -d "'" -f 2


for dir in */*.fits; do
    d=$( echo $dir | sed s"/.*VIMOS.//" | sed s"/T.*//" ) 
    t=$( imhead $dir | grep "HIERARCH ESO DPR TYPE" | cut -d "'" -f 2 )
    a=$( echo $dir | sed s"/\/.*//" )
    ob=${a: -1}
    if [ $t = "OBJECT" ]; then
        mv $dir vimos/$d/object$ob/
    fi
    if [ $t = "WAVE,LAMP" ]; then
        mv $dir vimos/$d/standard$ob/
    fi
    if [ $t = "FLAT,LAMP" ]; then
        mv $dir vimos/$d/standard$ob/
    fi
done


for dir in /Data/vimosindi/*/Bias/*/*[0-9][0-9][0-9].fits; do
    d=$( echo $dir | sed s"/.*VIMOS.//" | sed s"/T.*//" ) 
    a=$( echo $dir | sed s"/\/Bias.*//" )
    ob=${a: -1}
    cp $dir $d/bias/
done


for d in $( ls */VIMOS_IFU_OBS*.fits | cut -d"/" -f1 | sort -u ) ; do
    cd $d
    i=8
    for f in $( ls VIMOS_IFU_OBS*.fits ) ; do
        t=$( echo $i/8 | bc )
        mv $f $t/
        i=$(($i+1))
    done
    cd ..
done



for d in $( ls -d */standard[1-3]); do 
    cd $d
    Py3D vimos renameFiles 2012
    cd ../..
done

for d in $( ls -d */standard[1-3]); do 
    cd $d
    for f in $( ls VIMOS_IFU_LAMP*); do
        n=$( echo $f | cut -d"_" -f3 | sed 's/LAMP//' )
        o=$( echo $f | cut -d"_" -f4 )
        t=$( imhead $f | grep "DATE-OBS" | cut -d"T" -f3 | cut -d"'" -f1 )
        x=$( echo $t | awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3 }' )
        
        obs=$( ls ../../*/[1-3]/VIMOS_IFU_OBS$n* )
        tobs=$( imhead $obs |grep "DATE-OBS" | cut -d"T" -f3 | cut -d"'" -f1 )
        tobs=($tobs)
        xobs=$( echo $tobs | awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3 }' )
                        
        
        
    done
    cd ../..
done