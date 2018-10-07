#!/bin/bash

make

for i in {1..10}; do
   echo $i
   awk NR==$i a-CoO > tmp
#   awk NR==$i a1-NiO  > tmp
#   awk NR==$i a-FeO  > tmp
   cat tmp
   cat input tmp > seedsig.inp
   ./seedsig
   rm tmp
   mv sig-0102010100 sig-0102010100-$i
   mv sig-0202010100 sig-0202010100-$i
   mv sig-0302010100 sig-0302010100-$i
   mv sig-0402010100 sig-0402010100-$i
#   cp seedsig.log seedsig.log-$i
done 
