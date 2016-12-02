#!/bin/bash

enscf=$(grep "Final SCF" $1 | awk '{printf "%15.8f", $5}')
en1=$(grep 'E4(T)' $1 | awk '{s+=$3}END{print s}')
en2=$(grep 'E5(T)' $1 | awk '{s+=$3}END{print s}')
echo $enscf $(echo $en1 $en2 | awk '{printf "%15.8f\n", $1+$2}')

