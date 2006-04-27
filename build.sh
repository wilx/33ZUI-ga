#!/bin/sh

set -x 

gcc -O3 -W -Wall -std=c99 -pedantic -o ga ga.c -lm
