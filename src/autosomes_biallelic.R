# Filter for autosomes variants ----
autosomes <- alldata %>%
  filter(sample %in% qcpassing_samples) %>%
  filter(str_detect(chr, "chrX|chrM", negate = TRUE)) %>%
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
autosomes_mean_dp_variants <- alldata %>%
  filter(str_detect(chr, "chrX|chrM", negate = TRUE)) %>%
  filter(sample %in% qcpassing_samples) %>%
  group_by(variant) %>%
  summarize(
    n_individuals = n(),
    dp_mean = mean(dp, na.rm = TRUE),
    dp_sd = sd(dp, na.rm = TRUE),
    min_dp = min(dp, na.rm = TRUE),
    max_dp = max(dp, na.rm = TRUE)
  ) %>%
  mutate(var_qc_dp = ifelse(dp_mean >= 20, "pass", "fail"))

# Pivot autosomes dataframe ----
autosomes_pivoted_bi <- autosomes %>%
  group_by(hpos_surdez, variant, gt) %>%
  count() %>%
  pivot_wider(
    id_cols = c("hpos_surdez", "variant"),
    names_from = gt,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(call_rate = 1 - (`./.` / sum(across(where(is.numeric))))) %>%
  mutate(qc_call_rate = ifelse(call_rate >= 0.95, "pass", "fail")) %>%
  left_join({autosomes_mean_dp_variants %>%
      select(variant, var_qc_dp)}, by = "variant") %>%
  left_join({gene_names_alleles %>% select(variant, biallelic)}, by = "variant") %>%
  filter(biallelic == "yes") %>%
  rowwise() %>%
  mutate(hwep = hw_test_autosomes(`0/0`, `0/1`, `1/1`)) %>%
  ungroup()
