# *linc-mipep* and *linc-wrb* encode micropeptides that regulate chromatin accessibility in vertebrate-specific neural cells

- This repository contains R scripts and codes for various analyses of scRNA-seq and scATAC-seq data in the following preprint by Tornini et al.:
<https://www.biorxiv.org/content/10.1101/2022.07.21.501032v1>

- The input data are from Chromium Single Cell Multiome ATAC + Gene Expression by 10x Genomics (e.g. filtered_feature_bc_matrix.h5):
<https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression>

- Data processing and our customized analysis methods are based on or developed upon the Seurat toolkit:
<https://satijalab.org/seurat/index.html>

- The scripts and codes are shared without optimization in this repository, which may include analyses and results that are not reported in the preprint above. Repetition of scripts and codes was intended.

## File description
1. `README.md`
  : This current page
2. `sessionInfo.txt`
  : R session information from sessionInfo() in our analysis platform
3. `analysis.01.wt.R`
  : R scripts for pre-processing and QC of scRNA-seq and scATAC-seq data in WT
4. `analysis.02.mut.R`
  : R scripts for pre-processing and QC of scRNA-seq and scATAC-seq data in mutant
5. `analysis.03.merged.R`
  : R scripts and codes for various merged analyses of WT and mutant scRNA-seq and scATAC-seq data
6. `LICENSE`
  : MIT license

## Contacts
- Ho-Joon Lee, Ph.D.: **ho-joon.lee__at__yale.edu**
- Valerie Tornini, Ph.D.: **valerie.tornini__at__yale.edu**
- Antonio Giraldez, Ph.D.: **antonio.giraldez__at__yale.edu**

## License
Released under the MIT license. See LICENSE.
