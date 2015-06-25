#!/bin/bash

make
mpirun -n 8 --hostfile hostfile  ./ptsp
