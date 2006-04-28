#!/bin/sh

set -x 

gcc -O3 -save-temps -fverbose-asm -W -Wall -std=c99 -pedantic -DNDEBUG -o ga ga.c -lm
