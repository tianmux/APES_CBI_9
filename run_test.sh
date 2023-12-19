#!/bin/bash
make clean
make test
LOG_FILE="log_$(date +%Y-%m-%d_%H-%M-%S)_gcc.txt"
./build/testInputData ./inputs/input.json 0 > "$LOG_FILE" 2>&1
vim "$LOG_FILE" 

