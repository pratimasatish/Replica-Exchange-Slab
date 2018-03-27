#!/bin/bash --login

biaslist=( -0.750 -0.650 -0.600 -0.575 -0.550 -0.525 -0.500 -0.450 -0.375 -0.275 -0.200 -0.100 )

for bias in "${biaslist[@]}"
do
# 	echo ${bias}
	python less_frame.py -step 50 -bias ${bias} > sparse.${bias}
	./strip_dump.sh sparse.${bias} stripped.${bias}
	python colour_phase.py -f stripped.${bias} -div 0.77777778 > anti-gauche${bias}.xyz
# 	python colour_phase.py -f stripped.${bias} -div 0.77777778 > anti-gauche-old${bias}.xyz
	python colour_phase.py -f stripped.${bias} -div 0.575 > ord-disord${bias}.xyz
# 	python colour_phase.py -f stripped.${bias} -div 0.575 > ord-disord-old${bias}.xyz
        rm sparse.${bias} stripped.${bias}
done
