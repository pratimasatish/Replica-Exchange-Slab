#!/bin/bash --login

templist=( 345 348 350 351 352 355 )

for temp in "${templist[@]}"
do
# 	echo ${temp}
	python less_frame.py -step 50 -temp ${temp} > sparse.${temp}
	./strip_dump.sh sparse.${temp} stripped.${temp}
	python colour_phase.py -f stripped.${temp} -div 0.77777778 > anti-gauche-${temp}.xyz	
done
