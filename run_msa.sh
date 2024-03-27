#!/bin/sh
#echo input text    : $1
#echo input weights : $2
#echo utility U     : $3
#echo min len L     : $4
#echo create temporary files
./msa_create $1 $2 temp_$1_$2_$3_$4_text temp_$1_$2_$3_$4_weights
#echo run msa
/usr/bin/time -v ./msa temp_$1_$2_$3_$4_text temp_$1_$2_$3_$4_weights $3 $4
#echo delete temprary files temp_$1_$2_$3_$4_text temp_$1_$2_$3_$4_weights
rm temp_$1_$2_$3_$4_text temp_$1_$2_$3_$4_weights
#echo done
