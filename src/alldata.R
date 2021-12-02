# Gene names ----
gene_names <- read.xlsx(here("data_tables", "gene_names.xlsx"))

variants <- gene_names %>%
  left_join(
        {read.xlsx(here("data_tables", "variants.xlsx")) %>% mutate(chr = paste0("chr", Chromosome))}, 
        by = c("position" = "PositionVCF", "chr" = "chr")) %>%
  filter(Assembly == "GRCh38") %>%
  select(chr, position, gene_names, variant, Name, `RS#.(dbSNP)`, highestclinsig, ReferenceAlleleVCF, AlternateAlleleVCF) %>%
  distinct()

clinvar_pos <- read_delim(here("data_tables", "clinvar_pos.csv"),
                          col_names = c(
                            "chr",
                            "position",
                            "ref",
                            "alt"
                          )) %>%
  mutate(variant = glue("{chr}:{position}"))

clinvar_pos_multi_allelic <- clinvar_pos %>%
  group_by(variant) %>%
  count()

mt_vars <- gene_names %>% filter(chr == "chrMT") %>%
  select(-ref)

bi_vars <- gene_names %>% inner_join({clinvar_pos_multi_allelic %>% filter(n == 1)}, by = "variant") %>%
  filter(chr != "chrMT") %>%
  select(-c("n", "ref", "alt_list"))

multi_vars <- gene_names %>% inner_join({clinvar_pos_multi_allelic %>% filter(n > 1)}, by = "variant") %>%
  filter(chr != "chrMT") %>%
  select(-c("n", "ref", "alt_list"))

multi_vars_alleles <- clinvar_pos %>% 
  filter(variant %in% multi_vars$variant) %>%
  group_by(variant) %>%
  summarise(ref_list = paste0(ref, collapse = ","), alt_list = paste0(alt, collapse = ",")) %>%
  mutate(alt_alleles_n = str_count(alt_list, ",") + 1)

gene_names_alleles <- bind_rows(bi_vars, multi_vars) %>%
  left_join(multi_vars_alleles, by = "variant") %>%
  bind_rows(mt_vars) %>%
  mutate(biallelic = ifelse(str_detect(alt_list, ","), "no", "yes")) %>%
  mutate(chr_n = case_when(
    chr == "chrX" ~ 23,
    chr == "chrMT" ~ 25,
    chr != "chrX" | chr != "chrMT" ~ as.numeric(str_remove(chr, "chr"))
  )) %>%
  relocate(chr_n, .after = chr) %>%
  arrange(chr_n, position) %>%
  mutate(variant = str_replace(variant, "chrMT", "chrM")) %>%
  mutate(biallelic = replace_na(biallelic, "yes"))

# Sample metrics ----
metrics <- read.xlsx(here("data_tables", "sample_metrics.xlsx")) %>%
  mutate(s = str_remove(s, "-merge|-merg")) %>%
  filter(str_detect(routine, "_GRAR_")) %>%
  separate(routine, into = c(NA, "routine"), sep = "_GRAR_")

# Filenames ----
filenames <- read.xlsx(here("data_tables", "filenames_ids.xlsx")) %>%
  mutate(filenames2 = filenames) %>%
  separate(filenames2,into = c("sample", "routine"), sep = "_GRAR_") %>%
  mutate(sample = str_remove(sample, "-merge|-merg")) %>%
  mutate(routine = str_remove(routine, ".csv")) %>%
  filter(include == "yes")

load(here("data_tables", "probands.RData"))

hpos <- probands %>% 
  select(patient_id, phenotips_hpos) %>% 
  mutate(phenotips_hpos = str_remove_all(phenotips_hpos, '\\[|\\]|"')) %>%
  separate_rows(phenotips_hpos, sep = "}, ")%>%
  separate(phenotips_hpos, sep = ", ", into = c("hpo_id", "type", "label", "observed")) %>%
  mutate(hpo_id = str_remove(hpo_id, "\\{id: ")) %>%
  mutate(type = str_remove(type, "type: ")) %>%
  mutate(label = str_remove(label, "label: ")) %>%
  mutate(observed = str_remove(observed, "observed: "))

hpos_surdez <- hpos %>% filter(observed == "yes") %>%
  filter(str_detect(label, regex("deafness|hearing|deaf", ignore_case = TRUE))) %>%
  pull(patient_id) %>% unique()

filenames <- filenames %>% mutate(hpos_surdez = ifelse(patient_id %in% hpos_surdez, "yes", "no"))

# Sampleq QC ----
filenames_metrics <- filenames %>% left_join(metrics, by = c("sample" = "s", "routine" = "routine")) %>%
  mutate(contamination_score = ifelse(contamination_pct <= 2, 1, 0)) %>%
  mutate(q30bases_score = ifelse(q30bases >= 80*10^9, 1, 0)) %>%
  mutate(q30bases_pct_score = ifelse(q30bases_pct > 90, 1, 0)) %>%
  mutate(autos_coverage_score = ifelse(autosome_median_coverage >= 20, 1, 0)) %>%
  mutate(autos_call_pct_score = ifelse(autosome_callability_pct > 95, 1, 0)) %>%
  mutate(mapreads_pct_score = ifelse(mapped_reads_pct > 98, 1, 0)) %>%
  mutate(chimera_score = ifelse(chimeras_pct < 5, 1, 0)) %>%
  mutate(duplication_score = ifelse(duplication_pct < 10, 1, 0)) %>%
  mutate(insertsize_score = ifelse(median_insert_size > 300, 1, 0)) %>%
  mutate(score = rowSums(select(., ends_with("score")), na.rm = TRUE)) %>%
  mutate(qc = ifelse(score >= 6, "pass", "fail")) %>%
  mutate(qc = case_when(
    contamination_score == 0 ~ "fail",
    autos_coverage_score == 0 ~ "fail",
    score >= 6 ~ "pass",
    score < 6 ~ "fail"
  ))

qcpassing_samples <- filenames_metrics %>% filter(qc == "pass") %>%
  pull(sample)

# Read CSV files with genotypes ---
alldata <- lapply(seq_along(filenames_metrics$filenames), function(file) {
  
  read_delim(here("data_processed", filenames_metrics$filenames[file]), 
                      delim = " ",
                      col_names = c(
                        "chr",
                        "position",
                        "ref",
                        "alt",
                        "dp",
                        "gt"
                      ),
             col_types = cols(
               chr = "c",
               position = "i",
               ref = "c",
               alt = "c",
               dp = "i",
               gt = "c"
             )) %>%
    mutate(tmp = str_remove(filenames$filenames[file], ".csv")) %>%
    separate(tmp,into = c("sample", "routine"), sep = "_GRAR_") %>%
    mutate(sample = str_remove(sample, "-merge|-merg")) %>%
    mutate(variant = as.character(glue("{chr}:{position}"))) %>%
    relocate(sample) %>%
    relocate(variant, .after = sample)
  
}) %>% 
# Gather all data into a single dataframe
  bind_rows() %>%
# Fix genotype column
  mutate(gt = str_replace(gt, "\\|", "/"))
