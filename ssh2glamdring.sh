#!/bin/bash

# Bash script to ssh all necessary files to glamdring in one go.
reduced=/Data/vimos/cubes/$1.cube.combined.fits
analysis=/Data/vimos/analysis/$1
analysis1=$analysis/templates.txt
analysis2=$analysis/voronoi_2d_binning_output_$2.txt
analysis3=$analysis/voronoi_2d_binning_output2_$2.txt
p_file=/home/HOME/VIMOS_project/analysis/params.txt
g_file=/Data/vimos/analysis/galaxies.txt


#scp $reduced warrenj@glamdring.physics.ox.ac.uk:cubes/
scp -r $analysis1 $analysis2 $analysis3 warrenj@glamdring.physics.ox.ac.uk:analysis/$1/
scp $p_file warrenj@glamdring.physics.ox.ac.uk:
scp $g_file warrenj@glamdring.physics.ox.ac.uk:analysis/