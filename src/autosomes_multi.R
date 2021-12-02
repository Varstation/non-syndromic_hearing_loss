surdez_sample_size <- table(filenames_metrics$qc, filenames_metrics$hpos_surdez)["pass", "yes"]
nao_surdez_sample_size <- table(filenames_metrics$qc, filenames_metrics$hpos_surdez)["pass", "no"]

surdez_sex <- filenames_metrics %>% filter(include == "yes") %>% filter(qc == "pass")

n_females_surdez <- table(surdez_sex$sex, surdez_sex$hpos_surdez)["F","yes"]
n_males_surdez <- table(surdez_sex$sex, surdez_sex$hpos_surdez)["M","yes"]

n_females_nao_surdez <- table(surdez_sex$sex, surdez_sex$hpos_surdez)["F","no"]
n_males_nao_surdez <- table(surdez_sex$sex, surdez_sex$hpos_surdez)["M","no"]

# Filter for multiallele autosomes variants ----
sample_size <- length(unique(autosomes$sample))


# Filter for multiallele autosomes variants ----
autosomes_multi <- alldata %>%
  filter(sample %in% qcpassing_samples) %>%
  filter(str_detect(chr, "chrX|chrM", negate = TRUE)) %>%
  filter(variant %in% multi_vars_alleles$variant) %>%
  left_join({filenames %>% select(sample, hpos_surdez)}, by = "sample")

tmp <- autosomes_multi %>%
  group_by(hpos_surdez,variant, ref, alt, gt) %>%
  count() %>%
  pivot_wider(
    id_cols = c("hpos_surdez","variant", "ref", "alt"),
    names_from = gt,
    values_from = n,
    values_fill = 0
  ) %>%
  left_join({variants %>% select(-gene_names)}, by = c("variant" = "variant", "ref" = "ReferenceAlleleVCF", "alt" = "AlternateAlleleVCF")) %>%
  left_join({autosomes_mean_dp_variants %>%
      select(variant, var_qc_dp)}, by = "variant")

missing <- autosomes_multi %>%
  filter(alt == ".") %>%
  filter(gt == "./.") %>%
  distinct() %>%
  group_by(hpos_surdez, variant, gt) %>%
  count() %>%
  rename("./._new" = "n") %>%
  ungroup() %>%
  select(-gt)

homref <- autosomes_multi %>%
  filter(gt == "0/0") %>%
  distinct() %>%
  group_by(hpos_surdez, variant, gt) %>%
  count() %>%
  rename("0/0_new" = "n") %>%
  ungroup() %>%
  select(-gt)

corrected_genotypes_counts <- missing %>% 
  right_join({autosomes_multi %>% select(hpos_surdez, variant) %>% distinct()}, by = c("hpos_surdez" = "hpos_surdez", "variant" = "variant")) %>%
  right_join(homref, by = c("hpos_surdez" = "hpos_surdez", "variant" = "variant")) %>%
  mutate(`./._new` = replace_na(`./._new`, 0))

autosomes_pivoted_multi <- tmp %>%  
  left_join(corrected_genotypes_counts, by = c("hpos_surdez" = "hpos_surdez", "variant" = "variant")) %>%
  mutate(
    qc_call_rate = ifelse(
      hpos_surdez == "yes",
      
      ifelse( (1 - (`./._new` / surdez_sample_size)) >= 0.95 , "pass", "fail"),
      
      ifelse( (1 - (`./._new` / nao_surdez_sample_size)) >= 0.95 , "pass", "fail")
      
    )) %>%

  mutate(
    hwep = ifelse(
      hpos_surdez == "yes",
      ifelse(`1/2` == 0, hw_test_autosomes((surdez_sample_size - `0/1` - `1/1` - `./._new`), `0/1`, `1/1`),  hw_test_triallelic_autosome((surdez_sample_size - `0/1` - `1/1` - `1/2` - `./._new`), `0/1`, `1/1`, `1/2`)),
      ifelse(`1/2` == 0, hw_test_autosomes((nao_surdez_sample_size - `0/1` - `1/1` - `./._new`), `0/1`, `1/1`),  hw_test_triallelic_autosome((nao_surdez_sample_size - `0/1` - `1/1` - `1/2` - `./._new`), `0/1`, `1/1`, `1/2`))
    )) %>%
  
  distinct()

