#!/bin/bash --login

biaslist=( -0.750 -0.700 -0.675 -0.625 -0.600 -0.550 -0.500 -0.450 -0.425 -0.400 -0.375 -0.350 -0.275 -0.225 -0.175 -0.125 )
# biaslist=( -0.750 -0.700 -0.675 -0.625 -0.600 -0.550 -0.500 -0.450 -0.425 )
# biaslist=( -0.400 -0.375 -0.350 -0.275 -0.225 -0.175 -0.125 )

# for bias in "${biaslist[@]}"
# do
#     NLINES=$(echo $((5000*10322)) | bc)
# #     tail -n $NLINES dump.new.${bias}.xyz > dump.smaller.${bias}.xyz
#     split -l $NLINES dump.new.${bias}.xyz
#     mv xab dump.new.${bias}.xyz
#     rm xaa
# done

# for bias in "${biaslist[@]}"
# do
#     mv dump.smaller.${bias}.xyz dump.new.${bias}.xyz
# done

for bias in "${biaslist[@]}"
do
    NLINES=$(echo $((5000*10322)) | bc)
#     tail -n $NLINES dump.new.${bias}.xyz > dump.smaller.${bias}.xyz
    split -l $NLINES dump.${bias}.xyz
    mv xab dump.${bias}.xyz
    rm xaa
done

