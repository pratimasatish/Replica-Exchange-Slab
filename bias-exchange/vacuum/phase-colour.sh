#!/bin/bash --login

biaslist=( ord disord )

for bias in "${biaslist[@]}"
do
# 	echo ${bias}
	python less_frame.py -step 50 -bias ${bias} > sparse.${bias}
	./strip_dump.sh sparse.${bias} stripped.${bias}
	python colour_phase.py -f stripped.${bias} -div 0.77777778 > anti-gauche-${bias}.xyz
	python colour_phase.py -f stripped.${bias} -div 0.575 > ord-disord-${bias}.xyz
done
