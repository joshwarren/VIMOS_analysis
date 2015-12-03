#!/bin/bash

# Bash script to create gifs from kinematic maps produced by plot_results.py

for galaxy in ngc3557 ic1459 ic1531 ic4296 ngc0612 ngc1399 ngc3100 ngc7075\
    pks0718-34 eso443-g024
do
#    mkdir /Data/vimosindi/analysis/$galaxy/results/4200-/plots/\
#notinterpolated/gifs
    for plots in v sigma h3 h4
    do 
        files=(/Data/vimosindi/analysis/$galaxy/results/4200-/plots/\
notinterpolated/$plots*)

        s=$(identify ${files[*]} | awk '{print $3}' | sort -nr | head -n1)

        convert -delay 120 ${files[0]} -extent $s -delay 75 ${files[1]} -extent $s \
/Data/vimosindi/analysis/$galaxy/results/4200-/plots/notinterpolated/\
gifs/$plots.gif
    done
done