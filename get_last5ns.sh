#!/bin/bash --login

templist=( 345 348 350 351 352 355 )

for temp in "${templist[@]}"
do
	tail -n 110590000 dump.new.${temp}.xyz > dump.end.${temp}.xyz
done
