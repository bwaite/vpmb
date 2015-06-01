#!/bin/bash


for FILE  in $(find -name "vpm_decompression_input.json" | sort)
 do
    echo $FILE
    ../vpmb_c $FILE
    echo $?
done