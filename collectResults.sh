#!/bin/bash

#Bash script to find and collect the kinematic maps from the output of plot_results.py into a folder on the desktop.

for Type in v sigma h3 h4
do
    mkdir -p ~/Desktop/results/$Type-fields/gifs/

    for galaxy in ngc3557 ic1459 ic1531 ic4296 ngc0612 ngc1399 ngc3100 ngc7075\
        pks0718-34 eso443-g024
    do
        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/${Type}_field_4200-.png ~/Desktop/results/$Type-fields/$galaxy-$Type.png
        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/${Type}_uncert_field_4200-.png ~/Desktop/results/$Type-fields/$galaxy-${Type}_uncert.png


        cp /Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/gifs/${Type}.gif ~/Desktop/results/$Type-fields/gifs/$galaxy-${Type}.gif

    done
done