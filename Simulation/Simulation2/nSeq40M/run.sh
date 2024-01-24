#!/bin/bash

## Create slurm scripts for running the simulation
Rscript Create_slurm_jobs.R
sleep 60

## Submit the slurm jobs:
for FILENAME in Slurm_scripts/*
do
  sbatch $FILENAME
  sleep 900
done  

