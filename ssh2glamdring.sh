#!/bin/bash

# Bash script to ssh all necessary files to glamdring in one go.
reduced=/Data/vimos/cubes/$1.cube.combined.corr.fits
analysis=/Data/vimos/analysis/$1
analysis1=$analysis/templates.txt
analysis2=$analysis/$2/setup/voronoi_2d_binning_output.txt
analysis3=$analysis/$2/setup/voronoi_2d_binning_output2.txt
p_file=/home/HOME/VIMOS_project/analysis/params.txt
g_file=/Data/vimos/analysis/galaxies.txt

ssh warrenj@glamdring.physics.ox.ac.uk mkdir -p analysis/$1/$2/setup/
# scp $reduced warrenj@glamdring.physics.ox.ac.uk:cubes/
scp $analysis1 warrenj@glamdring.physics.ox.ac.uk:analysis/$1/
scp $analysis2 $analysis3 warrenj@glamdring.physics.ox.ac.uk:analysis/$1/$2/setup/
scp $p_file warrenj@glamdring.physics.ox.ac.uk:
scp $g_file warrenj@glamdring.physics.ox.ac.uk:analysis/
