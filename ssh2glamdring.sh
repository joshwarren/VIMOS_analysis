#!/bin/bash

# Bash script to ssh all necessary files to glamdring in one go.
reduced=/Data/vimosindi/reduced/$1
analysis=/Data/vimosindi/analysis/$1
analysis1=$analysis/templates.txt
analysis2=$analysis/voronoi_2d_binning_output.txt
analysis3=$analysis/voronoi_2d_binning_output2.txt
p_file=/home/warrenj/VIMOS_project/analysis/params.txt


scp -r $reduced warrenj@glamdring.physics.ox.ac.uk:reduced/
scp $analysis1 $analysis2 $analysis3 warrenj@glamdring.physics.ox.ac.uk:analysis/$1/
idl -e "create_params_file, '$1'"
scp $p_file warrenj@glamdring.physics.ox.ac.uk: