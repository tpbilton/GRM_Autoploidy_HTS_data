#!/bin/bash

## Create the slurm jobs for the simulations
Rscript Create_slurm_jobs.R
sleep 10

## ploidy 1
sbatch Slurm_scripts/Sim_poly_GUSrelate_p1_d5_e0.sl
sleep 3000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p1_d5_e1.sl
sleep 3000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p1_d5_e10.sl
sleep 3000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p1_d25_e0.sl
sleep 3000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p1_d25_e1.sl
sleep 3000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p1_d25_e10.sl
sleep 3000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p1_d50_e0.sl
sleep 3000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p1_d50_e1.sl
sleep 3000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p1_d50_e10.sl
sleep 3000

# ploidy 4
sbatch Slurm_scripts/Sim_poly_GUSrelate_p4_d5_e0.sl
sleep 12000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p4_d5_e1.sl
sleep 12000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p4_d5_e10.sl
sleep 12000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p4_d25_e0.sl
sleep 12000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p4_d25_e1.sl
sleep 12000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p4_d25_e10.sl
sleep 12000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p4_d50_e0.sl
sleep 12000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p4_d50_e1.sl
sleep 12000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p4_d50_e10.sl
sleep 12000

## Ploidy 2
sbatch Slurm_scripts/Sim_poly_GUSrelate_p2_d5_e0.sl
sleep 6000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p2_d5_e1.sl
sleep 6000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p2_d5_e10.sl
sleep 6000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p2_d25_e0.sl
sleep 6000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p2_d25_e1.sl
sleep 6000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p2_d25_e10.sl
sleep 6000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p2_d50_e0.sl
sleep 6000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p2_d50_e1.sl
sleep 6000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p2_d50_e10.sl
sleep 6000

## Ploidy 3
sbatch Slurm_scripts/Sim_poly_GUSrelate_p3_d5_e0.sl
sleep 9000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p3_d5_e1.sl
sleep 9000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p3_d5_e10.sl
sleep 9000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p3_d25_e0.sl
sleep 9000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p3_d25_e1.sl
sleep 9000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p3_d25_e10.sl
sleep 9000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p3_d50_e0.sl
sleep 9000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p3_d50_e1.sl
sleep 9000
sbatch Slurm_scripts/Sim_poly_GUSrelate_p3_d50_e10.sl
sleep 9000



