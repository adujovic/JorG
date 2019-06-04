#!/bin/bash

i=0
files=""
for d in $@; do
    touch tmp.$i
    echo $d > tmp.$i
    grep -A 88 "magnetization (x)" $d | grep "\." | awk '{a[$1]=$5} END{for(i=1; i<=NR; i++) if(i in a) print a[i]}' >> tmp.$i
    files+=" tmp.$i"
    (( i+=1 ))
done

echo $files
paste $files
rm $files
exit
