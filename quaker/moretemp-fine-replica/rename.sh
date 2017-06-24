#!/bin/bash --login

biaslist=$(seq -0.8500 0.025 0.1000)

for bias in ${biaslist}
do
	mv theta${bias}.txt theta-355.${bias}.txt
	mv theta_x${bias}.txt theta_x-355.${bias}.txt
done
