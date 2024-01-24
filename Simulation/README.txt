Folder contains R scripts for performing the two simulations in the manuscript "Construction of relatedness matrices in autopolyploid populations using low depth high-throughput sequencing data".

- Simulation 1: Performance of the different GRMs when the mean read depth and sequencing error is varied for different ploidy levels. To run this simulation:
  1. sh run.sh
  2. Rscript Sim_poly_analysis.R #only run this after all the slurm jobs from one have finished.

- Simulation 2:Performance of GUSrelate for a fixed sequencing effort.
  There are two folders:
  * nSeq10M: Simulation using a total sequencing effort of 10M reads
  * nSeq40M: Simulation using a total sequencing effort of 40M reads
  The structure of the files in each folder is the same. Within each folder, the following commands need to be excuted in order:
  1. $ sh run.sh
  2. $ Rscript Sim_poly_analysis.R #only run this after all the slurm jobs from one have finished.
  Once both of these commands have finished in each folder, then go into the first level of the Simulation2 folder and run the command:
  * $ Rscript Plot_Figure4.R

##### Notes: ######
To run pedigreeSim, the following files/folder from the PedigreeSim folder have to be copied into the Simulation1 and Simulation2 folders:
- File: PedigreeSim.jar
- Folder: lib