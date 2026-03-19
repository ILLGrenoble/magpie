#!/bin/bash

mcstas --verbose -o test.c test.instr
gcc -I.. -L../build -o test test.c -lmagpie_c
