#!/bin/bash --login

biaslist=( -0.750 -0.700 -0.675 -0.625 -0.600 -0.550 -0.500 -0.450 -0.425 -0.400 -0.375 -0.350 -0.275 -0.225 -0.175 -0.125 )

for bias in "${biaslist[@]}"
do
# 	python zcos_th.py -bias ${bias} -show -reshape 
	python zcos_th.py -bias ${bias} -reshape 
done

# for bias in ${biaslist}
# do
# 	python xcos_th.py -bias ${bias}
# done
