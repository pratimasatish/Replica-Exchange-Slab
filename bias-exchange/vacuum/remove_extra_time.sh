#!/bin/bash --login

NLINES=$(echo $((10000*10322)) | bc)
split -l $NLINES dump-disord-moremoremore-400.xyz
mv xab dump-moremoremore-last${bias}.xyz
mv xaa dump-moremoremore-first${bias}.xyz

