#!/bin/bash --login

templist=$(seq 355 3 358)

for temp in ${templist}
do
# 	echo ${temp}
	awk '(/^2/ || /^1/) && $3 >= 0 && $3 <= 50' dump.new.${temp}.xyz > lig.${temp}
done
