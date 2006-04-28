#!/bin/sh

set -x

echo 
for i in `seq 1 100` ; do
    ./ga -l -g 50000 -f |tee -a f101-ordered.txt
done
cp f101-ordered.txt f101-ordered-`date +%Y%m%d%H%M%S`.txt

for i in `seq 1 100` ; do
    ./ga -s -l -g 50000 -f |tee -a f101-shuffled.txt
done
cp f101-shuffled.txt f101-shuffled-`date +%Y%m%d%H%M%S`.txt
