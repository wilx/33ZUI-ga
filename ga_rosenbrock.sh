#!/bin/sh

set -x

echo 
for i in `seq 1 100` ; do
    ./ga -l -g 50000 -f |tee -a ros-ordered.txt
done
cp ros-ordered.txt ros-ordered-`date +%Y%m%d%H%M%S`.txt

for i in `seq 1 100` ; do
    ./ga -s -l -g 50000 -f |tee -a ros-shuffled.txt
done
cp ros-shuffled.txt ros-shuffled-`date +%Y%m%d%H%M%S`.txt
