#!/bin/bash --login

templist=( 345 348 350 351 352 355 )

for temp in "${templist[@]}"
do
	python zcos_th.py -temp ${temp} -show
done
