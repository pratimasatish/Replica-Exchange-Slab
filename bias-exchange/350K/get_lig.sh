#!/bin/bash --login

biaslist=( -0.750 -0.650 -0.600 -0.575 -0.550 -0.525 -0.500 -0.450 -0.375 -0.275 -0.200 -0.100 )

for bias in "${biaslist[@]}"
do
# 	awk '(/^2/ || /^1/) && $3 >= 0' dump.${bias}.xyz > lig.old.${bias}
	awk '(/^2/ || /^1/) && $3 >= 0' dump.new.${bias}.xyz > lig.${bias}
done
