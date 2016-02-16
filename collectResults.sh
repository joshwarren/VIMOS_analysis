#!/bin/bash

#Bash script to find and collect the kinematic maps from the output of plot_results.py into a folder on the desktop.

for galaxy in $(find /Data/vimosindi/analysis/ -mindepth 1 -maxdepth 1 -type d | sed 's/\/Data\/vimosindi\/analysis\///')
do
    for Type in vel sigma h3 h4
    do
    stellar_dir=~/Desktop/results/stellar/$Type-fields/
    gas_dir=~/Desktop/results/hot_gas/$Type-fields/
    comparison_dir=~/Desktop/results/comparison/$Type-fields

    mkdir -p $stellar_dir/gifs/
    mkdir -p $gas_dir/gifs/
    mkdir -p $comparison_dir/


        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/stellar_${Type}_field_4200-.png $stellar_dir/$galaxy-$Type.png
        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/stellar_${Type}_uncert_field_4200-.png $stellar_dir/$galaxy-${Type}_uncert.png


        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/gifs/stellar_${Type}.gif $stellar_dir/gifs/$galaxy-${Type}.gif


        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/gas_${Type}_field_4200-.png $gas_dir/$galaxy-$Type.png
        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/gas_${Type}_uncert_field_4200-.png $gas_dir/$galaxy-${Type}_uncert.png


        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/gifs/gas_${Type}.gif $gas_dir/gifs/$galaxy-${Type}.gif



        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/gifs/stellar-gas_${Type}.gif $comparison_dir/$galaxy-${Type}.gif


    done
    mkdir -p ~/Desktop/results/residuals/

    cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/median_residual_4200-.png ~/Desktop/results/residuals/$galaxy-residuals.png


    mkdir -p ~/Desktop/results/chi2/

    cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/chi2_4200-.png ~/Desktop/results/chi2/$galaxy-Chi2.png

done

