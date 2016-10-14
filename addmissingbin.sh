#!/bin/bash
gal=eso443-g024
bin='454' # existing bin
bin2='455' # missing bin
#whichtype=pop_MC # for stellar populations
whichtype=gas_MC # for kinematics
base_dir=''

# if doesn't exist, copy the wavelength file first (then is exempt from setting all values to 0)
if [ ! -f $base_dir/Data/vimos/analysis/$gal/$whichtype/lambda/$bin2.dat ]
then
	echo $base_dir/Data/vimos/analysis/$gal/$whichtype/lambda/$bin2.dat
	cp $base_dir/Data/vimos/analysis/$gal/$whichtype/lambda/$bin.dat $base_dir/Data/vimos/analysis/$gal/$whichtype/lambda/$bin2.dat
fi

# Copies template number correctly, but give zero weighting
if [ ! -f $base_dir/Data/vimos/analysis/$gal/$whichtype/temp_weights/$bin2.dat ]
then 
    echo $base_dir/Data/vimos/analysis/$gal/$whichtype/temp_weights/$bin2.dat 
    sed s'/ [0-0].[0-9]*/ 0.0/g' $base_dir/Data/vimos/analysis/$gal/$whichtype/temp_weights/$bin.dat > temp
    mv temp $base_dir/Data/vimos/analysis/$gal/$whichtype/temp_weights/$bin2.dat 
fi



for f in $( find $base_dir/Data/vimos/analysis/$gal/$whichtype/ -print | grep /$bin.dat )
do
    f2=$( echo $f | sed s"/$bin/$bin2/" )
    # If file does not already exist
    if [ ! -f $f2 ]
    then
    	echo $f2
    	cp $f $f2
    	sed s'/[0-9]/0/g' $f2 > temp 
    	mv temp $f2
    	# Undo changing OIII name in templates file.
    	# NB: stellar template numbers will be set to 0, but have 0 weighting anyway
    	sed s'/\[OIII\]0000d/\[OIII\]5007d/' $f2 > temp
    	mv temp $f2
    fi
done