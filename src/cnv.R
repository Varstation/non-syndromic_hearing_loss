# Read ANNOTSV column names vector ----
source(here("src", "annotsv_colnames.R"))

# Fisher exact test helper function ----
# source(here("functions", "fisher_wrap.R"))

# Read CNV files ----
cnv_files <- list.files(here("cnv_annot"), recursive = TRUE, pattern = "_cnv_cnv_annotted.txt.tsv")

cnv_df <- data.frame(filename = cnv_files) %>%
  mutate(filename2 = filename) %>%
  separate(filename2, into = c("routine", "sample"), sep = "/") %>%
  mutate(routine = str_remove(routine, "GRAR_")) %>%
  mutate(sample = str_remove(sample, "_cnv_cnv_annotted.txt.tsv")) %>%
  mutate(sample = str_remove(sample, "-merge|-merg")) %>%
  left_join(filenames_metrics, by = c("sample" = "sample", "routine" = "routine")) %>%
  filter(include == "yes", qc == "pass")

allcnv <- lapply(seq_along(cnv_df$filename), function(file) {
  
  routine <- str_remove(base::strsplit(cnv_df$filename[file], split = "/")[[1]][1], "GRAR_")
  
  sample_name <- str_remove_all(base::strsplit(cnv_df$filename[file], split = "/")[[1]][2], ".txt|-merge|-merg|_cnv_cnv_annotted.txt.tsv")
  
  tmp <- readr::read_tsv(here("cnv_annot", cnv_df$filename[file]),
                         col_types = "cciidcccccicccccccicicciciiciidciidccccddcccdcccccccddcccccccccccccddcccccdddddddcccccic") %>%
    mutate(
      sample = sample_name,
      routine = routine)
  
}) %>% bind_rows()

# All CNVs (`AnnotSV type` == "full") -----
allcnv_fulls <- allcnv %>% 
  filter(`AnnotSV type` == "full") %>% 
  left_join({cnv_df %>% select(sample, routine, hpos_surdez, sex, status)}, 
            by = c("sample" = "sample", "routine" = "routine"))  %>%
  separate(Einstein, into = c("genotype", "sm", "cn", "bc", "pe"), sep = ":") %>%
  mutate(cn = as.numeric(cn)) %>%
  distinct() %>% 
  mutate(sex_check = case_when(
    `SV chrom` == "X" & sex == "M" & cn == 1 ~ "remove",
    `SV chrom` == "X" & sex == "F" & cn == 2 ~ "remove",
    `SV chrom` == "X" & sex == "M" & cn != 1 ~ "keep",
    `SV chrom` == "X" & sex == "F" & cn != 2 ~ "keep",
    `SV chrom` != "X" ~ "keep"
  )) %>%
  filter(sex_check == "keep") %>%
  mutate(genes = ifelse(str_detect(`Gene name`, "/"), glue("multiple genes (chr{`SV chrom`})"), `Gene name`)) %>%
  mutate(hard_filter = ifelse(FILTER == "PASS", "pass", "fail")) %>%
  filter(`AnnotSV ranking` >= 4) %>%
  relocate(sample)

## Multiloci deletions ----
multiloci_deletions <- allcnv_fulls %>% filter(str_detect(`Gene name`, "/")) %>% filter(ALT =="<DEL>") %>%
  relocate(sample, routine, `Gene name`, location) %>%
  group_by(hpos_surdez,`SV chrom`,`Gene name`, `SV length`) %>% count() %>%
  arrange(desc(n)) 

multiloci_deletions_rename <- multiloci_deletions %>%
  rename("Hearing loss status" = "hpos_surdez", 
         "Chromosome" = "SV chrom",
         "Loci" = "Gene name",
         "Deletion length (bp)" = "SV length", 
         "Individuals affected" = "n")

## Multiloci duplications ----
multiloci_duplications <- allcnv_fulls %>% filter(str_detect(`Gene name`, "/")) %>% filter(ALT =="<DUP>") %>%
  relocate(sample, routine, `Gene name`, location) %>%
  group_by(hpos_surdez,`SV chrom`,`Gene name`, `SV length`) %>% count() %>%
  arrange(desc(n)) 

multiloci_duplications_rename <- multiloci_duplications %>%
  ungroup() %>%
  select(-`SV length`) %>%
  rename("Hearing loss status" = "hpos_surdez", 
         "Chromosome" = "SV chrom",
         "Loci" = "Gene name",
         "Individuals affected" = "n")
  
# All CNVs (`AnnotSV type` == "split") -----
allcnv_splits <- allcnv %>% 
  filter(`AnnotSV type` == "split") %>% 
  left_join({cnv_df %>% select(sample, routine, hpos_surdez, sex, status)}, 
            by = c("sample" = "sample", "routine" = "routine"))  %>%
  separate(Einstein, into = c("genotype", "sm", "cn", "bc", "pe"), sep = ":") %>%
  mutate(cn = as.numeric(cn)) %>%
  distinct() %>% 
  mutate(sex_check = case_when(
    `SV chrom` == "X" & sex == "M" & cn == 1 ~ "remove",
    `SV chrom` == "X" & sex == "F" & cn == 2 ~ "remove",
    `SV chrom` == "X" & sex == "M" & cn != 1 ~ "keep",
    `SV chrom` == "X" & sex == "F" & cn != 2 ~ "keep",
    `SV chrom` != "X" ~ "keep"
  )) %>%
  filter(sex_check == "keep") %>%
  mutate(genes = ifelse(str_detect(`Gene name`, "/"), glue("multiple genes (chr{`SV chrom`})"), `Gene name`)) %>%
  mutate(hard_filter = ifelse(FILTER == "PASS", "pass", "fail")) %>%
  filter(`AnnotSV ranking` >= 4) %>%
  relocate(sample)

splits_genes_surdez <- allcnv_splits %>% 
  filter(`Gene name` %in% unique(gene_names$gene_names))

splits_genes_surdez_samples <- splits_genes_surdez %>% select(sample, hpos_surdez) %>% distinct()
table(splits_genes_surdez_samples$hpos_surdez)

## Location coordinates ----
location_coordinates <- splits_genes_surdez %>% select(ALT, `SV chrom`,`Gene name`, `SV start`, `SV end`, location) %>% distinct() %>%
  group_by(ALT, `SV chrom`,`Gene name`, location) %>%
  summarise(start = min(`SV start`), end = max(`SV end`))

## Sample size -----
surdez_sample_size <- table(filenames_metrics$qc, filenames_metrics$hpos_surdez)["pass", "yes"]
nao_surdez_sample_size <- table(filenames_metrics$qc, filenames_metrics$hpos_surdez)["pass", "no"]

surdez_sex <- filenames_metrics %>% filter(include == "yes") %>% filter(qc == "pass") %>% filter(!(sample %in% sex_gt_disagreement))

n_females_surdez <- table(surdez_sex$sex, surdez_sex$hpos_surdez)["F","yes"]
n_males_surdez <- table(surdez_sex$sex, surdez_sex$hpos_surdez)["M","yes"]

n_females_nao_surdez <- table(surdez_sex$sex, surdez_sex$hpos_surdez)["F","no"]
n_males_nao_surdez <- table(surdez_sex$sex, surdez_sex$hpos_surdez)["M","no"]

## Autosomal deletions grouped by location ----
alldels_autosomes_v2 <- splits_genes_surdez %>%
  filter(`SV chrom` != "X") %>%
  filter(ALT == "<DEL>") %>%
  group_by(hpos_surdez,
           ALT,
           `SV chrom`,
           `Gene name`,
           location, 
           cn, 
           `AnnotSV ranking`, 
           frameshift,
           Phenotypes) %>% count() %>%
  mutate(chr_n = as.numeric(str_replace(`SV chrom`, "X", "23"))) %>%
  arrange(chr_n) %>%
  select(-chr_n) %>%
  pivot_wider(!n, names_from = c("hpos_surdez", "cn"), values_from = n, values_fill = 0) %>%
  mutate(yes_2 = (surdez_sample_size - yes_0 - yes_1)) %>%
  mutate(no_2 = (nao_surdez_sample_size - no_0 - no_1)) %>%
  mutate(del_yes = (2 *(yes_0)) + yes_1 ) %>%
  mutate(del_no = (2 * (no_0)) + no_1) %>%
  mutate(wt_yes = (2 *(yes_2)) + yes_1 ) %>%
  mutate(wt_no = (2 * (no_2)) + no_1) %>%
  mutate(alt_yes_freq = del_yes / (2 * surdez_sample_size)) %>%
  mutate(alt_no_freq = del_no / (2 * nao_surdez_sample_size)) %>%
  mutate(hgf_no = no_1 / nao_surdez_sample_size) %>%
  mutate(hgf_yes = yes_1 / surdez_sample_size)

## X chromosome deletions grouped by location ----
alldels_x_v2 <- splits_genes_surdez %>%
  filter(`SV chrom` == "X") %>%
  filter(ALT == "<DEL>") %>%
group_by(sex,
    hpos_surdez,
           ALT,
           `SV chrom`,
           `Gene name`,
           location, 
           cn, 
           `AnnotSV ranking`, 
    frameshift,
           Phenotypes) %>%
  count() %>%
  pivot_wider(id_cols = c("sex", "ALT", `SV chrom`, `Gene name`, "location", `AnnotSV ranking`, frameshift, Phenotypes) , names_from = c("hpos_surdez", "cn"), values_from = n, values_fill = 0) %>%
  mutate(alt_no_freq = no_1 / ((2 * n_females_nao_surdez) + n_males_nao_surdez))  %>%
  mutate(alt_yes_freq = 0) %>%
  mutate(hgf_no = no_1 / n_females_nao_surdez) %>%
  mutate(hgf_yes = 0) %>%
  mutate(yes_0 = 0, yes_1 = 0, no_0 = 0)

### Publication Deletion table ----
deletions_table <- 
alldels_autosomes_v2 %>%
  mutate(sex = "any") %>%
  select(ALT, sex, `SV chrom`, `Gene name`, location, `AnnotSV ranking`, frameshift, Phenotypes, yes_0, yes_1, no_0, no_1, alt_yes_freq, hgf_yes, alt_no_freq, hgf_no) %>%
  bind_rows(
    {
      alldels_x_v2 %>% select(ALT, sex, `SV chrom`, `Gene name`, location, `AnnotSV ranking`, frameshift, Phenotypes,yes_0, yes_1, no_0, no_1, alt_yes_freq, hgf_yes, alt_no_freq, hgf_no)
    }
    
  ) %>%
  left_join(location_coordinates, by = c("ALT" = "ALT", "SV chrom" = "SV chrom", "Gene name" = "Gene name", "location" = "location")) %>%
  ungroup() %>%
  mutate(chrom_pos = glue("chr{`SV chrom`}:{start}-{end}")) %>%
  relocate(chrom_pos, .after = ALT) %>%
  select(-c(`SV chrom`, start, end))

deletions_table_rename <- deletions_table %>%
  rename(
    "CNV Type" = "ALT",
    "Sex" = "sex",
         "Chromosome coordinate ranges" = "chrom_pos",
         "Gene region" = "location",
         "Frameshift" = "frameshift",
         "Copy number = 0 (hearing loss group)" = "yes_0",
         "Copy number = 1 (hearing loss group)" = "yes_1",
         "Copy number = 0 (no hearing loss group)" = "no_0",
         "Copy number = 1 (no hearing loss group)" = "no_1",
         "Alternate allele frequency (hearing loss group)" = "alt_yes_freq",
         "Heterozygote genotypic frequency (hearing loss group)" = "hgf_yes",
         "Alternate allele frequency (no hearing loss group)" = "alt_no_freq",
         "Heterozygote genotypic frequency (no hearing loss group)" = "hgf_no"
         )
 
## Autosomal duplications (grouped by location) ----
alldups_autosomes_v2 <- splits_genes_surdez %>%
  filter(`SV chrom` != "X") %>%
  filter(ALT == "<DUP>") %>%
  group_by(hpos_surdez,
           ALT,
           `SV chrom`,
           `Gene name`,
           location, 
           cn, 
           `AnnotSV ranking`, 
           frameshift,
           Phenotypes) %>% count() %>%
  mutate(chr_n = as.numeric(str_replace(`SV chrom`, "X", "23"))) %>%
  arrange(chr_n) %>%
  select(-chr_n) %>%
  pivot_wider(!n, names_from = c("hpos_surdez", "cn"), values_from = n, values_fill = 0) %>%
  ungroup() %>%
  # mutate(no_2 = nao_surdez_sample_size - rowSums(select(., starts_with("no_")))) %>%
  # mutate(yes_2 = surdez_sample_size - rowSums(select(., starts_with("yes_")))) %>%
  mutate(alt_yes_freq = (yes_3) / (2 * surdez_sample_size)) %>%
  mutate(alt_no_freq = (no_3) / (2 * nao_surdez_sample_size)) %>%
  mutate(hgf_yes = (yes_3) /  surdez_sample_size) %>%
  mutate(hgf_no = (no_3) / nao_surdez_sample_size)

duplications_table_autosomes <- alldups_autosomes_v2 %>%
  mutate(sex = "any") %>%
  select(ALT, sex, `SV chrom`, `Gene name`, location, `AnnotSV ranking`, frameshift, Phenotypes, yes_3, no_3, no_4, no_5, alt_yes_freq, hgf_yes, alt_no_freq, hgf_no) %>%
  left_join(location_coordinates, by = c("ALT" = "ALT", "SV chrom" = "SV chrom", "Gene name" = "Gene name", "location" = "location")) %>%
  ungroup() %>%
  mutate(chrom_pos = glue("chr{`SV chrom`}:{start}-{end}")) %>%
  relocate(chrom_pos, .after = ALT) %>%
  select(-c(`SV chrom`, start, end))

### Publication autosomal Duplication table ---- 
duplications_table_autosomes_rename <- duplications_table_autosomes %>%
  rename(
    "CNV Type" = "ALT",
    "Sex" = "sex",
    "Chromosome coordinate ranges" = "chrom_pos",
    "Gene region" = "location",
    "Frameshift" = "frameshift",
    "Copy number = 3 (hearing loss group)" = "yes_3",
    "Copy number = 3 (no hearing loss group)" = "no_3",
    "Copy number = 4 (no hearing loss group)" = "no_4",
    "Copy number = 5 (no hearing loss group)" = "no_5",
    "Alternate allele frequency (hearing loss group)" = "alt_yes_freq",
    "Heterozygote genotypic frequency (hearing loss group)" = "hgf_yes",
    "Alternate allele frequency (no hearing loss group)" = "alt_no_freq",
    "Heterozygote genotypic frequency (no hearing loss group)" = "hgf_no"
  )

## X chromosome duplications (grouped by location) ----

alldups_x_v2 <- splits_genes_surdez %>%
  filter(`SV chrom` == "X") %>%
  filter(ALT == "<DUP>") %>%
  group_by(hpos_surdez,
           sex,
           ALT,
           `SV chrom`,
           `Gene name`,
           location, 
           cn, 
           `AnnotSV ranking`, 
           frameshift,
           Phenotypes) %>%
  relocate(sex, `Gene name`, location, cn, GD_AF, Phenotypes, .after = sample) %>%
  count() %>%
  pivot_wider(id_cols = c("sex", "ALT", `SV chrom`, `Gene name`, "location", `AnnotSV ranking`, frameshift, Phenotypes) , names_from = c("hpos_surdez", "cn"), values_from = n, values_fill = 0) %>%
  mutate(extra_cn_no = no_3 + no_2) %>%
  mutate(extra_cn_yes = 0)

### Publication Duplication in X chromosome table ----
xdups_allele_freq <- alldups_x_v2 %>% group_by(`Gene name`) %>% summarise(affected_x = sum(extra_cn_no)) %>% mutate(alt_freq_yes = 0, alt_freq_no = affected_x / ((2 * n_females_nao_surdez) + n_males_nao_surdez))

xdups_final <- alldups_x_v2 %>% left_join(xdups_allele_freq, by = "Gene name") %>%
  select(-c("no_3", "no_2", "affected_x")) %>%
  relocate(extra_cn_yes, .before = extra_cn_no) %>%
  left_join(location_coordinates, by = c("ALT" = "ALT", "SV chrom" = "SV chrom", "Gene name" = "Gene name", "location" = "location")) %>%
  ungroup() %>%
  mutate(chrom_pos = glue("chr{`SV chrom`}:{start}-{end}")) %>%
  relocate(chrom_pos, .after = ALT) %>%
  select(-c(`SV chrom`, start, end))

xdups_final_rename <- xdups_final %>%
  rename(
    "Sex" = "sex",
    "CNV Type" = "ALT",
    "Chromosome coordinate ranges" = "chrom_pos",
    "Gene region" = "location",
    "Frameshift" = "frameshift",
    "Individuals with extra copies in X chromosome (hearing loss group)" = "extra_cn_yes",
    "Individuals with extra copies in X chromosome (no hearing loss group)" = "extra_cn_no",
    "Duplication alleles frequency (hearing loss group)" = "alt_freq_yes",
    "Duplication alleles frequency (no hearing loss group)" = "alt_freq_no"
  )


# Export publication tables ----

wb <- createWorkbook()

addWorksheet(wb, sheetName = "Multiloci DELS", gridLines = TRUE)

writeData(wb, sheet = "Multiloci DELS", x = multiloci_deletions_rename)

addWorksheet(wb, sheetName = "Multiloci DUPS", gridLines = TRUE)

writeData(wb, sheet = "Multiloci DUPS", x = multiloci_duplications_rename)

addWorksheet(wb, sheetName = "Deletions", gridLines = TRUE)

writeData(wb, sheet = "Deletions", x = deletions_table_rename)

addWorksheet(wb, sheetName = "Duplications (autos.)", gridLines = TRUE)

writeData(wb, sheet = "Duplications (autos.)", x = duplications_table_autosomes_rename)

addWorksheet(wb, sheetName = "Duplications (chr X)", gridLines = TRUE)

writeData(wb, sheet = "Duplications (chr X)", x = xdups_final_rename)

saveWorkbook(wb, here("output_tables", "CNV results.xlsx"), overwrite = TRUE)


 
