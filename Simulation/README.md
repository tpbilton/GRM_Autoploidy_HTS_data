# Instructions for Running the Simulations

Simulation files associated with the manuscript by Bilton et al. (2024).

## Folder structure:

-   Simulation 1: Performance of the different GRMs when the mean read depth and sequencing error is varied for different ploidy levels. To run this simulation:

    1.  Run the command:
```
$ sh run.sh
```

    2.  When all the code from step 1 has finished (e.g., all the slurm jobs have finished), run the command:
```
$ Rscript Sim_poly_analysis.R
```

-   Simulation 2: Performance of GUSrelate for a fixed sequencing effort. There are two folders:

    -   nSeq10M: Simulation using a total sequencing effort of 10M reads
    -   nSeq40M: Simulation using a total sequencing effort of 40M reads The structure of the files in each folder is the same. Within each folder, the following commands need to be excuted in order:

    1. Run the command:  
```
$ sh run.sh
```

    2. When all the slurm jobs from step one have finished, run the command:
```
$ Rscript Sim_poly_analysis.R
``` 

    Once step 1 and 2 are complete in both the nSeq10M and nSeq40M folders, then go into the root level of the Simulation2 folder and run the command:
```
$ Rscript Plot_Figure4.R
```
- PedigreeSim: Copy of the the code for the software PedigreeSim used in the simulations.

#### R package versions:

The version of GUSbase and GUSrelate that were used in the simulations can be downloaded using the following code:

```
    devtools::install_github("tpbilton/GUSbase", ref = "92119b9")
    devtools::install_github("tpbilton/GUSrelate", ref = "d93411f")
```

Note: the `ref` argument indicates which commit from GitHub to install the packages from.

A full list of the R packages and version numbers used for the simulation are shown below:

```
    R version 4.0.3 (2020-10-10)
    Platform: x86_64-conda-linux-gnu (64-bit)
    Running under: CentOS Linux 7 (Core)

    Matrix products: default
    BLAS/LAPACK: /bifo/scratch/Potato_GBS_Teagasc/conda/lib/libopenblasp-r0.3.12.so

    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
     [9] LC_ADDRESS=C               LC_TELEPHONE=C
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

    attached base packages:
    [1] parallel  stats     graphics  grDevices utils     datasets  methods
    [8] base

    other attached packages:
    [1] doParallel_1.0.16 iterators_1.0.13  foreach_1.5.1     AGHmatrix_2.1.3
    [5] GUSrelate_0.2.0   plotly_4.9.2.2    ggplot2_3.3.3     GUSbase_0.2.1

    loaded via a namespace (and not attached):
     [1] pillar_1.4.7      compiler_4.0.3    tools_4.0.3       digest_0.6.27
     [5] lattice_0.20-41   jsonlite_1.7.2    lifecycle_0.2.0   tibble_3.0.4
     [9] gtable_0.3.0      viridisLite_0.3.0 pkgconfig_2.0.3   rlang_0.4.10
    [13] Matrix_1.3-2      withr_2.3.0       dplyr_1.0.2       httr_1.4.2
    [17] generics_0.1.0    vctrs_0.3.6       htmlwidgets_1.5.3 gbRd_0.4-11
    [21] grid_4.0.3        tidyselect_1.1.0  glue_1.4.2        data.table_1.13.6
    [25] R6_2.5.0          Rdpack_2.1        purrr_0.3.4       tidyr_1.1.2
    [29] magrittr_2.0.1    scales_1.1.1      codetools_0.2-18  ellipsis_0.3.1
    [33] htmltools_0.5.0   rbibutils_2.0     colorspace_2.0-0  lazyeval_0.2.2
    [37] munsell_0.5.0     crayon_1.3.4      zoo_1.8-8
```

#### Java version:

The simulations were performed using the following Java version

```
    openjdk version "1.8.0_282"
    OpenJDK Runtime Environment (build 1.8.0_282-b08)
    OpenJDK 64-Bit Server VM (build 25.282-b08, mixed mode)
```

**Note**
To run pedigreeSim, the following files/folder from the PedigreeSim folder have to be copied into the Simulation1 and Simulation2 folders:

- File: PedigreeSim.jar
- Folder: lib

### Funding

-   Ministry of Business, Innovation and Employment via its funding of the "Genomics for Production & Security in a Biological Economy" programme (Contract ID C10X1306).
-   Ministry of Business, Innovation and Employment via its Strategic Science Investment Fund (SSIF) to AgResearch.

### References

Bilton, T.P., Sharma, S.K., Schofield, M.R., Black, M.A., Jacobs, J.M.E., Bryan, G.J. & Dodds, K.G. (2024). Construction of relatedness matrices in autopolyploid populations using low depth high-throughput sequencing data. Unpublished Manuscript
