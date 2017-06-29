#!/bin/bash --login

templist=$(seq $1 $2 $3)

for temp in ${templist}
do
	python zcos_th.py -temp ${temp}
done

for temp in ${templist}
do
	python xcos_th.py -temp ${temp}
done
