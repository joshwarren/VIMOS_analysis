#grep r1 multirun.sh.o48019 > glamdring_result.txt

#cat glamdring_result.txt | awk -f ~/VIMOS_project/analysis/awkFile.txt > glamdring_result5.txt

#grep ^r1 glamdring_result5.txt > glamdring_result6.txt

#cat glamdring_result6.txt | sed 's/r[1-2] *$//' | sed 's/[a-z]+ *$//' > glamdring_result7.txt


grep r1 multirun.sh.o55487 | awk -f ~/VIMOS_project/analysis/awkFile.txt | grep ^r1 | sed 's/r[1-2] *$//' | sed 's/[a-z]* *$//'| sed 's/[A-Z]*$//' > glamdring_result.txt

grep r2 multirun.sh.o55487 | awk -f ~/VIMOS_project/analysis/awkFile2.txt | grep ^r2 | sed 's/r[1-2] *$//' | sed 's/[a-z]* *$//' | sed 's/[A-Z]*$//' > glamdring_result2.txt
