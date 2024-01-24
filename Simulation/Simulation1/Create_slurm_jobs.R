### Script for creating slurm jobs for running the simulation
## If additional slurm setting are required, then added them in the file below:

slurm_dir = "Slurm_scripts/"
dir.create(slurm_dir)

## set up the different parameters for the simulation study
ploidLevels = c(1,2,3,4)
Depths = c(5,25,50)
epsilon = c(0,0.001,0.01)

## create the slurm file and run
for(ploid in ploidLevels){
  for(depth in Depths){
    for(ep in epsilon){
      slurm_file = paste0(slurm_dir,"Sim_poly_GUSrelate_p",ploid,"_d",depth,"_e",ep*1000,".sl")
      cat("#!/bin/bash -e",
          paste0("#SBATCH --job-name=poly_sim_p",ploid,"_d",depth,"_e",ep*1000),
          "#SBATCH --ntasks=1",
          "#SBATCH --cpus-per-task=1",
          "#SBATCH --mem-per-cpu=1GB",
          "#SBATCH --time 7-24",
          paste0("#SBATCH --chdir=", getwd()),
          "#SBATCH --output=SlurmOut/Slurm_%x_%j.out",
          "#SBATCH --error=SlurmOut/Slurm_%x_%j.err",
          "#SBATCH --array=1-500%50",
          "",
          paste("srun Rscript Sim_poly_GUSrelate.R $SLURM_ARRAY_TASK_ID", ploid,depth,ep),
          file = slurm_file, sep="\n")
    }
  }
}
  
