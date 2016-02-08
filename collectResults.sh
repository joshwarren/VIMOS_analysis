#!/bin/bash

#Bash script to find and collect the kinematic maps from the output of plot_results.py into a folder on the desktop.

for Type in vel sigma h3 h4
do
    stellar_dir=~/Desktop/results/stellar/$Type-fields/
    gas_dir=~/Desktop/results/hot_gas/$Type-fields/
    comparison_dir=~/Desktop/results/comparison/$Type-fields

    mkdir -p $stellar_dir/gifs/
    mkdir -p $gas_dir/gifs/
    mkdir -p $comparison_dir/



    for galaxy in ngc3557 ic1459 ic1531 ic4296 ngc0612 ngc1399 ngc3100 ngc7075\
        pks0718-34 eso443-g024
    do
        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/stellar_${Type}_field_4200-.png $stellar_dir/$galaxy-$Type.png
        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/stellar_${Type}_uncert_field_4200-.png $stellar_dir/$galaxy-${Type}_uncert.png


        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/gifs/stellar_${Type}.gif $stellar_dir/gifs/$galaxy-${Type}.gif


        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/gas_${Type}_field_4200-.png $gas_dir/$galaxy-$Type.png
        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/gas_${Type}_uncert_field_4200-.png $gas_dir/$galaxy-${Type}_uncert.png


        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/gifs/gas_${Type}.gif $gas_dir/gifs/$galaxy-${Type}.gif



        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/gifs/stellar-gas_${Type}.gif $comparison_dir/$galaxy-${Type}.gif
    done
done



mkdir -p ~/Desktop/results/residuals/
for galaxy in ngc3557 ic1459 ic1531 ic4296 ngc0612 ngc1399 ngc3100 ngc7075\
    pks0718-34 eso443-g024
do
    cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/median_residual_4200-.png ~/Desktop/results/residuals/$galaxy-residuals.png

done




