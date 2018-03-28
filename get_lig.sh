#!/bin/bash --login

templist=( 345 348 350 351 352 355 )

for temp in "${templist[@]}"
do
# 	echo ${temp}
	awk '(/^2/ || /^1/) && $3 >= 0 && $3 <= 50' dump.end.${temp}.xyz > lig.${temp}
done
