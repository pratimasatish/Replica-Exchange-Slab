#!/bin/bash --login

templist=$(seq 355 3 358)

for temp in ${templist}
do
	python zcos_th.py -temp ${temp} -show
done
