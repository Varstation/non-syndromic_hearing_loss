library(here)
library(openxlsx)
library(tidyverse)
library(glue)
library(htmltab)

# Source HWE test function ----
source(here("functions", "hw_test.R"))

# 1- Gather data and perform sample QC ----
source(here("src","alldata.R"))

# 2- Variant QC and data tabulation and ----

## 2.1 Autosome variants ----
### 2.1.1 Biallelic variants ----
source(here("src","autosomes_biallelic.R"))

### 2.1.2 Multi-allele variants ----
source(here("src","autosomes_multi.R"))

## 2.2 Chr X variants ----
source(here("src","chrx.R"))

## 2.3 Chr M variants ----
source(here("src", "chrm.R"))

# 3- Export ClinVar variants tables
source(here("src", "to_excel.R"))

# 4- CNV
source(here("src", "cnv.R"))

# 5- SV
source(here("src", "sv.R"))
