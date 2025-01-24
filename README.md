# MetabolSSMF

Simplex-structure matrix factorisation: application of soft clustering to metabolomics data

The package is built upon the 'Biobase' package from 'Bioconductor.' Before installing 'MetabolSSMF,' please ensure that 'Biobase' is installed by following the provided instructions.

# Install the package

```{r}
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("Biobase")

devtools::install_github("WenxuanLiu1996/MetabolSSMF")
```

# Load the package

```{r}
library(MetabolSSMF)
```

