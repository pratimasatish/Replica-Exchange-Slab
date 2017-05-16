#!/bin/bash --login

templist=$(seq $1 $2 $3)

for temp in ${templist}
do
	python xcos_th.py -temp ${temp} -start 0 -show
done
