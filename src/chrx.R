# Source HWE test function ----
chrx <- alldata %>%
  filter(sample %in% qcpassing_samples) %>%
  filter(str_detect(chr, "chrX")) %>%
  pivot_wider(
    id_cols = sample,
    names_from = variant,
    values_from = gt) %>%
  mutate(
    across(everything(), ~replace_na(.x, "./."))
  ) %>%
  pivot_longer(!sample, names_to = "variant", values_to = "gt") %>%
  left_join({filenames %>% select(sample, sex)}, by = "sample") %>%
  left_join({filenames %>% select(sample, hpos_surdez)}, by = "sample")

# Sample QC: sex/genotype disagreement ----
hemiz_vars_females <- chrx %>% filter(gt == "0" | gt == "1", sex == "F") %>% 
  pull(sample) %>% unique()

homoz_vars_males <- chrx %>% filter(gt == "0/0", sex == "M") %>% 
  pull(sample) %>% unique()

sex_gt_disagreement <- c(hemiz_vars_females, homoz_vars_males)

# Variant QC ----
## variant coverage ----
chrx_mean_dp_variants <- 
  alldata %>%
  filter(sample %in% qcpassing_samples) %>%
  filter(str_detect(chr, "chrX")) %>%
  filter(!sample %in% sex_gt_disagreement) %>%
  left_join({filenames %>% select(sample, sex)}, by = "sample") %>%
  group_by(variant) %>%
  summarize(
    n_individuals = n(),
    dp_mean = mean(dp, na.rm = TRUE),
    dp_sd = sd(dp, na.rm = TRUE),
    min_dp = min(dp, na.rm = TRUE),
    max_dp = max(dp, na.rm = TRUE)
  ) %>%
  mutate(var_qc_dp = ifelse(dp_mean >= 20, "pass", "fail"))

# Pivot chrx dataframe ----
chrx_pivoted <- chrx %>%
  filter(!sample %in% sex_gt_disagreement) %>%
  group_by(hpos_surdez, sex, variant, gt) %>%
  count() %>%
  pivot_wider(
    id_cols = c("hpos_surdez", "variant"),
    names_from = c("sex", "gt"),
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(call_rate_f = 1 - (`F_./.` / (`F_0/0` + `F_0/1` + `F_./.`))) %>%
  mutate(heterozygote_freq = (`F_0/1`) / (`F_0/0` + `F_0/1`)) %>%
  mutate(hemizygote_freq = (`M_1`) / (`M_0` + `M_1`)) %>%
  mutate(call_rate_m = 1) %>%
  mutate(ref_freq = ((2 * `F_0/0`) + `F_0/1` + M_0) / (2 * (`F_0/0` + `F_0/1`) + M_0)) %>%
  mutate(alt1_freq = 1 - ref_freq) %>%
  mutate(qc_call_rate = ifelse(call_rate_f >= 0.95, "pass", "fail")) %>%
  left_join({chrx_mean_dp_variants %>%
      select(variant, var_qc_dp)}, by = "variant") %>%
  left_join({gene_names_alleles %>% select(variant, biallelic)}, by = "variant") %>%
  mutate(`F_1/1` = 0) %>%
  relocate(`F_1/1`, .after = `F_0/1`) %>%
  relocate(`F_./.`, .before = `F_0/0`) %>%
  rowwise() %>%
  mutate(hwep = hw_test_xlinked(M_0, M_1, `F_0/0`, `F_0/1`, `F_1/1`)) %>%
  ungroup()
