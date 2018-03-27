#!/bin/bash --login

biaslist=( -0.750 -0.700 -0.675 -0.625 -0.600 -0.550 -0.500 -0.450 -0.425 -0.400 -0.375 -0.350 -0.275 -0.225 -0.175 -0.125 )

for bias in "${biaslist[@]}"
do
# 	echo ${bias}
	python less_frame.py -step 50 -bias ${bias} > sparse.${bias}
	./strip_dump.sh sparse.${bias} stripped.${bias}
# 	python colour_phase.py -f stripped.${bias} -div 0.77777778 > anti-gauche${bias}.xyz
# 	python colour_phase.py -f stripped.${bias} -div 0.77777778 > anti-gauche-old${bias}.xyz
# 	python colour_phase.py -f stripped.${bias} -div 0.575 > ord-disord${bias}.xyz
	python colour_phase.py -f stripped.${bias} -div 0.575 > ord-disord-old${bias}.xyz
done
