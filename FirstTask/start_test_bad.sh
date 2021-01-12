#!/bin/bash

sudo docker run -v $PWD:/sandbox:Z abouteiller/mpi-ft-ulfm mpirun --oversubscribe -mca btl tcp,self -np 10 a.out 1
