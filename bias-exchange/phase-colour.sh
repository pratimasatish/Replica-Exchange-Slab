#!/bin/bash --login

biaslist=( 360 375 )

for bias in "${biaslist[@]}"
do
	python less_frame.py -step 50 -bias ${bias} > sparse.${bias}
# 	awk '(/^2/ || /^1/) && $3 >= -1' sparse.${bias} > sparselig.${bias}
# 	python calc-ang-per-time.py -bias ${bias}
	./strip_dump.sh sparse.${bias} stripped.${bias}
	python colour_phase.py -f stripped.${bias} -div 0.77777778 > anti-gauche${bias}.xyz
 	python colour_phase.py -f stripped.${bias} -div 0.575 > ord-disord${bias}.xyz
	python colour_phase.py -f stripped.${bias} -div 0.25 > DO-div-${bias}.xyz
done
