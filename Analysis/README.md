# Instructions for Running the Analysis

Script required to reproduce the analysis of the Potato dataset in the manuscript by Bilton et al. (2024).

There is only one one script that needs to be run: **Pot_analysis_comp.R**.

Note that the following supplementary files needed to be download and put in the same folder as the script:

- Supplementary_File3.csv
- Supplementary_file4.tab.ra
- Supplementary_File5.csv
- Supplementary_file6.tab.ra

### R versions

The version of GUSbase and GUSrelate that were used in the simulations can be downloaded using the following code:

```
    devtools::install_github("tpbilton/GUSbase", ref = "02d570d")
    devtools::install_github("tpbilton/GUSrelate", ref = "0898175")
```

R session information when running the analysis:

```
R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=English_New Zealand.utf8  LC_CTYPE=English_New Zealand.utf8    LC_MONETARY=English_New Zealand.utf8
[4] LC_NUMERIC=C                         LC_TIME=English_New Zealand.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] MethComp_1.30.0 AGHmatrix_2.1.4 GUSrelate_0.2.2 plotly_4.10.4   ggplot2_3.4.4   GUSbase_0.2.3  

loaded via a namespace (and not attached):
  [1] colorspace_2.1-0        ellipsis_0.3.2          class_7.3-22            flextable_0.9.4         fs_1.6.3               
  [6] httpcode_0.3.0          rstudioapi_0.15.0       proxy_0.4-27            remotes_2.4.2.1         fansi_1.0.6            
 [11] lubridate_1.9.3         xml2_1.3.6              codetools_0.2-19        splines_4.2.1           cachem_1.0.8           
 [16] knitr_1.45              pkgload_1.3.3           jsonlite_1.8.8          shiny_1.8.0             compiler_4.2.1         
 [21] httr_1.4.7              Matrix_1.6-4            fastmap_1.1.1           lazyeval_0.2.2          cli_3.6.2              
 [26] later_1.3.2             htmltools_0.5.7         tools_4.2.1             coda_0.19-4             gtable_0.3.4           
 [31] glue_1.6.2              dplyr_1.1.4             Rcpp_1.0.11             fontquiver_0.2.1        vctrs_0.6.5            
 [36] crul_1.4.0              nlme_3.1-164            iterators_1.0.14        xfun_0.41               stringr_1.5.1          
 [41] rbibutils_2.2.16        ps_1.7.5                timechange_0.2.0        mime_0.12               miniUI_0.1.1.1         
 [46] lifecycle_1.0.4         devtools_2.4.5          zoo_1.8-12              scales_1.3.0            ragg_1.2.7             
 [51] promises_1.2.1          fontLiberation_0.1.0    curl_5.2.0              memoise_2.0.1           pander_0.6.5           
 [56] gdtools_0.3.5           stringi_1.8.2           fontBitstreamVera_0.1.1 desc_1.4.3              foreach_1.5.2          
 [61] e1071_1.7-14            pkgbuild_1.4.3          zip_2.3.0               epiR_2.0.67             Rdpack_2.6             
 [66] rlang_1.1.2             pkgconfig_2.0.3         systemfonts_1.0.5       evaluate_0.23           lattice_0.22-5         
 [71] purrr_1.0.2             sf_1.0-14               htmlwidgets_1.6.4       processx_3.8.3          tidyselect_1.2.0       
 [76] magrittr_2.0.3          R6_2.5.1                generics_0.1.3          profvis_0.3.8           DBI_1.2.1              
 [81] pillar_1.9.0            withr_2.5.2             units_0.8-5             survival_3.5-7          tibble_3.2.1           
 [86] crayon_1.5.2            gfonts_0.2.0            uuid_1.2-0              KernSmooth_2.23-22      utf8_1.2.4             
 [91] rmarkdown_2.25          officer_0.6.3           urlchecker_1.0.1        usethis_2.2.2           grid_4.2.1             
 [96] data.table_1.14.10      callr_3.7.3             digest_0.6.33           classInt_0.4-10         xtable_1.8-4           
[101] tidyr_1.3.0             httpuv_1.6.13           textshaping_0.3.7       openssl_2.1.1           munsell_0.5.0          
[106] viridisLite_0.4.2       BiasedUrn_2.0.11        sessioninfo_1.2.2       askpass_1.2.0          
```

### Funding

-   Ministry of Business, Innovation and Employment via its funding of the "Genomics for Production & Security in a Biological Economy" programme (Contract ID C10X1306).
-   Ministry of Business, Innovation and Employment via its Strategic Science Investment Fund (SSIF) to AgResearch.

### References

Bilton, T.P., Sharma, S.K., Schofield, M.R., Black, M.A., Jacobs, J.M.E., Bryan, G.J. & Dodds, K.G. (2024). Construction of relatedness matrices in autopolyploid populations using low depth high-throughput sequencing data. Unpublished Manuscript
