#!/bin/bash

Rscript rstan_simu_globloclaplace.R SIMU_H1 14

HPS="SIMU_H0 SIMU_H1"
#DATASETS="1 97 98 99 100"
for d in $(seq 15 100);
do
    for h in $HPS
    do
        # Rscript rstan_simu_laplace.R $h $d
        Rscript rstan_simu_globloclaplace.R $h $d
    done
    done
