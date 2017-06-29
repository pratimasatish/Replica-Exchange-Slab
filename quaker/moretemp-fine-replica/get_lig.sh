#!/bin/bash --login

templist=$(seq $1 $2 $3)

for temp in ${templist}
do
	awk '(/^2/ || /^1/) && $3 >= 0' dump.new.${temp} > lig.${temp}
done
