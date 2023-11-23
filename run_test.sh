#!/bin/bash
#make clean
#make test
make
clear
./build/testInputData ./inputs/input.json 0 >out
vim out

