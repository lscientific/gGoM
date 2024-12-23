# Supplementary Codes to paper ``Generalized Grade-of-Membership Estimation for High-dimensional Locally Dependent Data''

## Utility files
- `gom.R' includes the proposed method 
- `gibbs_util.R' includes the Gibbs sampling utility functions for simulation I

## Simulation experiments
- `sim_flatten.R' corresponds to Simulation I and generates Table 1, Table 2, Figure 2
- `sim_loc.R' corresponds to Simulation II and generates Figure 3
- `sim_pois.R' corresponds to Simulation III and generates Figure 4


## Three real data applications
- `hapmap.R', `ANES.R', `single cell.R' correspond to the three real data analyses
- The pre-processed HapMap3 data can be downloaded from [https://figshare.com/s/9b4d5964af498d167e85]. and the raw data can be accessed from [https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3]
- The ANES 2022 pilot study data is downloaded from [https://electionstudies.org/data-center/]
- The single cell data can be accessed by running
```
install_github('kkdey/singleCellRNASeqMouseDeng2014') 
library(singleCellRNASeqMouseDeng2014)
```
