#!/bin/bash

make
mpirun -n 4 ./ptsp
