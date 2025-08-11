# R Analysis Container

This container provides R with required packages for the h3achipimputation pipeline's analysis and visualization processes.

## Included R Packages

- **Core tidyverse packages**: dplyr, tidyr, readr, stringr, purrr, tibble
- **Data processing**: data.table
- **Visualization**: ggplot2, scales, gridExtra, cowplot, viridis, RColorBrewer

## Usage

This container is automatically built by GitHub Actions and used by the following pipeline processes:
- `average_r2` - Calculate average R-squared values
- `plot_r2_SNPcount` - Plot imputed SNPs vs mean R2
- `plot_hist_r2_SNPcount` - Plot histograms of SNPs vs R2
- `plot_MAF_r2` - Plot MAF vs R2
- `filter_info_by_target` - Filter imputation info by target
- `plot_freq_comparison` - Plot frequency comparisons
- `plot_r2_SNPpos` - Plot R2 by SNP position

## Container Registry

After push to GitHub, the container will be available at:
```
ghcr.io/mamana/r-analysis:1.0
```

## Local Testing

To build locally for testing:
```bash
docker build -t mamana/r-analysis:1.0 .
```

To convert to Singularity:
```bash
singularity pull docker://ghcr.io/mamana/r-analysis:1.0
```