#!/bin/bash --login

biaslist=( -0.750 -0.700 -0.650 -0.625 -0.600 -0.580 -0.560 -0.540 -0.520 -0.500 -0.425 -0.350 -0.275 -0.225 -0.175 -0.125 )

for bias in "${biaslist[@]}"
do
# 	python zcos_th.py -bias ${bias} -reshape
	python zcos_th.py -bias ${bias} -reshape -show 
done

# for bias in ${biaslist}
# do
# 	python xcos_th.py -bias ${bias}
# done
