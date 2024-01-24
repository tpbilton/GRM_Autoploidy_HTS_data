### Script for creating slurm jobs for running the optimal depth simulation
## If additional slurm setting are required, then added them in the file below:

slurm_dir = "Slurm_scripts/"
dir.create(slurm_dir)

## set up the different parameters for the simulation study
ploidLevels = c(1,2,3,4)
Depths = c(1:30,35,40)

## create the slurm file and run
for(ploid in ploidLevels){
  for(depth in Depths){
    slurm_file = paste0(slurm_dir,"Sim_poly_depth_GUSrelate_p",ploid,"_d",depth,".sl")
    cat("#!/bin/bash -e",
        paste0("#SBATCH --job-name=poly_sim_opd_p",ploid,"_d",depth),
        "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=1",
        "#SBATCH --mem-per-cpu=4GB",
        "#SBATCH --time 7-24",
        paste0("#SBATCH --chdir=",getwd()),
        "#SBATCH --error=out/slurm-%A-%a.err",
        "#SBATCH --output=out/slurm-%A-%a.out",
        "#SBATCH --array=1-500%50",
        "",
        paste("srun Rscript Sim_poly_depth_GUSrelate.R $SLURM_ARRAY_TASK_ID", ploid, depth),
        file = slurm_file, sep="\n")
  }
}
  
