#!/bin/bash --login

templist=$(seq $1 $2 $3)

for temp in ${templist}
do
	python zcos_th.py -temp ${temp} -start 0 -show
done
