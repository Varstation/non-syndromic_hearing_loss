chrm <- alldata %>%
  filter(sample %in% qcpassing_samples) %>%
  filter(str_detect(chr, "chrM"))%>%
  pivot_wider(
    id_cols = sample,
    names_from = variant,
    values_from = gt) %>%
  mutate(
    across(everything(), ~replace_na(.x, "./."))
  ) %>%
  pivot_longer(!sample, names_to = "variant", values_to = "gt") %>%
  left_join({filenames %>% select(sample, hpos_surdez)}, by = "sample")

# Variant QC ----
## variant coverage ----
chrm_mean_dp_variants <- 
  alldata %>%
  filter(sample %in% qcpassing_samples) %>%
  filter(str_detect(chr, "chrM")) %>%
  group_by(variant) %>%
  summarize(
    n_individuals = n(),
    dp_mean = mean(dp, na.rm = TRUE),
    dp_sd = sd(dp, na.rm = TRUE),
    min_dp = min(dp, na.rm = TRUE),
    max_dp = max(dp, na.rm = TRUE)
  ) %>%
  mutate(var_qc_dp = ifelse(dp_mean >= 20, "pass", "fail"))

# Pivot chrm dataframe ----
chrm_pivoted <- chrm %>%
  group_by(hpos_surdez, variant, gt) %>%
  count() %>%
  pivot_wider(
    id_cols = c("hpos_surdez", "variant"),
    names_from = gt,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(call_rate = 1 - (`./.` / sum(across(where(is.numeric))))) %>%
  left_join({chrm_mean_dp_variants %>%
      select(variant, var_qc_dp)}, by = "variant") %>%
  left_join({gene_names_alleles %>% select(variant, biallelic, gene_names, chr_n, position)}, by = "variant") %>%
  relocate(gene_names) %>%
  arrange(position) %>%
  select(-position)