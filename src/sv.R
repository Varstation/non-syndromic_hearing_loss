# Read CNV files ----
sv_files <- list.files(here("sv_annot"), recursive = TRUE, pattern = "_sv_annotted.txt.tsv")

sv_df <- data.frame(filename = sv_files) %>%
  mutate(filename2 = filename) %>%
  separate(filename2, into = c("routine", "sample"), sep = "/") %>%
  mutate(routine = str_remove(routine, "GRAR_")) %>%
  mutate(sample = str_remove(sample, "_sv_annotted.txt.tsv")) %>%
  mutate(sample = str_remove(sample, "-merge|-merg")) %>%
  left_join(filenames_metrics, by = c("sample" = "sample", "routine" = "routine")) %>%
  filter(include == "yes", qc == "pass")

# All SV split/filtered ----
allsv_split_filtered <- lapply(seq_along(sv_df$filename), function(file) {
  
  routine <- str_remove(base::strsplit(sv_df$filename[file], split = "/")[[1]][1], "GRAR_")
  
  sample_name <- str_remove_all(base::strsplit(sv_df$filename[file], split = "/")[[1]][2], ".txt|-merge|-merg|_sv_annotted.txt.tsv")
  
  tmp <- readr::read_tsv(here("sv_annot", sv_df$filename[file]),
                         skip = 1, 
                         col_names = annotsv_colnames,
                         col_types = "cciidcccccicccccccicicciciiciidciidccccddcccdcccccccddcccccccccccccddcccccdddddddcccccic") %>%
    mutate(
      sample = sample_name,
      routine = routine) %>%
    filter(FILTER == "PASS") %>%
    filter(`AnnotSV type` == "split") %>%
    filter(`Gene name` %in% unique(gene_names$gene_names)) %>%
    filter(`AnnotSV ranking` >= 4) %>%
    select(sample, 
           routine, 
           `SV chrom`,
           `SV start`,
           `SV end`,
           `SV length`,
           `SV type`,
           `Gene name`,
           REF,
           ALT,
           Einstein,
           location, 
           location2,
           distNearestSS,
           nearestSStype,
           tx,
           `CDS length`,
           frameshift,
           `tx length`,
           GD_AF,
           Phenotypes,
           `AnnotSV ranking`,
           `ranking decision criteria`
           ) %>% 
    left_join({sv_df %>% select(sample, routine, hpos_surdez, sex, status)}, 
              by = c("sample" = "sample", "routine" = "routine")) %>%
    separate(Einstein, into = c("gt", "ft", "gq", "pl", "pr"), sep = ":") %>%
    distinct()

}) %>% bind_rows()

## Location coordinates ----
sv_unique_samples <- allsv_split_filtered %>% select(hpos_surdez, sample, `SV type`, `SV chrom`,`Gene name`, location, gt,`AnnotSV ranking`, frameshift, Phenotypes) %>% distinct()

sv_location_coordinates <- allsv_split_filtered %>% select(`SV type`, `SV chrom`,`Gene name`, `SV start`, `SV end`, location) %>% distinct() %>%
  group_by(`SV type`, `SV chrom`,`Gene name`, location) %>%
  summarise(start = min(`SV start`), end = max(`SV end`))

## Known variants ----

shearer_table1 <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4067994/table/T1/"

previous_cnvs <- htmltab(shearer_table1)

dels_ins_dups <- allsv_split_filtered %>% 
  filter(`Gene name` %in% previous_cnvs$Gene) %>% 
  filter(frameshift == "yes" | location == "txStart-txEnd") %>% 
  filter(`SV type` != "INV")

dels_ins_dups_autosomes <-
dels_ins_dups %>%
  filter(`SV chrom` != "X") %>%
  group_by(
    hpos_surdez,
    `SV type`,
    `SV chrom`,
    `Gene name`,
    location, 
    gt, 
    `AnnotSV ranking`, 
    frameshift,
    Phenotypes) %>%
  count() %>%
  mutate(chr_n = as.numeric(str_replace(`SV chrom`, "X", "23"))) %>%
  arrange(chr_n) %>%
  select(-chr_n) %>%
  pivot_wider(!n, names_from = c("hpos_surdez", "gt"), values_from = n, values_fill = 0) %>%
  relocate(`yes_0/1`, .before = `no_0/1`) %>%
  mutate(alt_yes_freq = `yes_0/1` / (2 * surdez_sample_size)) %>%
  mutate(alt_no_freq = `no_0/1` / (2 * nao_surdez_sample_size)) %>%
  mutate(hgf_yes = `yes_0/1` / surdez_sample_size) %>%
  mutate(hgf_no = `no_0/1` / nao_surdez_sample_size) %>%
  arrange(desc(alt_no_freq))

dels_ins_dups_x <-
  dels_ins_dups %>%
  filter(`SV chrom`== "X") %>%
  group_by(
    sex,
    hpos_surdez,
    `SV type`,
    `SV chrom`,
    `Gene name`,
    location, 
    gt, 
    `AnnotSV ranking`, 
    frameshift,
    Phenotypes) %>%
  count() %>%
pivot_wider(id_cols = c("sex", "SV type", `SV chrom`, `Gene name`, "location", `AnnotSV ranking`, frameshift, Phenotypes), names_from = c("hpos_surdez", "gt"), values_from = n, values_fill = 0) %>%
  mutate(alt_yes_freq = 0) %>%
mutate(alt_no_freq = `no_0/1` / ((2 * n_females_nao_surdez) + n_males_nao_surdez)) %>%
  mutate(hgf_yes = 0) %>%
  mutate(hgf_no = `no_0/1` / n_females_nao_surdez)


dels_ins_dups_final <- dels_ins_dups_autosomes %>%
  mutate(sex = "any") %>%
  relocate(sex) %>%
  bind_rows(dels_ins_dups_x) %>%
  select(`SV type`, `SV chrom`, sex, `Gene name`, location, 
         `AnnotSV ranking`, frameshift, Phenotypes, `yes_0/1`, `no_0/1`, alt_yes_freq, hgf_yes, alt_no_freq, hgf_no) %>%
  left_join(sv_location_coordinates, by = c("SV type" = "SV type", "SV chrom" = "SV chrom", "Gene name" = "Gene name", "location" = "location")) %>%
  ungroup() %>%
  mutate(chrom_pos = glue("chr{`SV chrom`}:{start}-{end}")) %>%
  relocate(chrom_pos, .after = `SV type`) %>%
  select(-c(`SV chrom`, start, end))

dels_ins_dups_final_rename <- dels_ins_dups_final %>%
  rename(
    "SV Type" = "SV type",
    "Chromosome coordinate ranges" = "chrom_pos",
    "Sex" = "sex",
    "Gene region" = "location",
    "Frameshift" = "frameshift",
    "N heterozygotes (0/1) (hearing loss group)" = "yes_0/1",
    "N heterozygotes (0/1) (no hearing loss group)" = "no_0/1",
    "Alternate allele frequency (hearing loss group)" = "alt_yes_freq",
    "Heterozygote genotypic frequency (hearing loss group)" = "hgf_yes",
    "Alternate allele frequency (no hearing loss group)" = "alt_no_freq",
    "Heterozygote genotypic frequency (no hearing loss group)" = "hgf_no"
  )

wb <- createWorkbook()

addWorksheet(wb, sheetName = "SV results", gridLines = TRUE)

writeData(wb, sheet = "SV results", x = dels_ins_dups_final_rename)

saveWorkbook(wb, here("output_tables", "SV results.xlsx"), overwrite = TRUE)
  