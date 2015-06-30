#!/bin/bash

pnum=$1
if [ -z $pnum ]; then
	pnum=16
fi
time mpirun -n $pnum --hostfile hostfile  ./ptsp
