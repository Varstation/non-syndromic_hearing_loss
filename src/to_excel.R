# Output annovar ----
annovar <- read_tsv(here("data_tables", "input_annovar.txt.hg38_multianno.txt"), na = c("N/A", ".")) %>%
  select(dbsnp, Chr, Start, Ref, Alt, CLNSIGCONF, gnomAD_genome_ALL, starts_with("ABRaOM"), REVEL_score, REVEL_rankscore, SIFT_score, SIFT_pred, SIFT_converted_rankscore, starts_with("Polyphen2"), PVS1, PS1, PS2, PS3, PS4, PM1, PM2, PM3, PM4, PM5, PM6, PP1, PP2, PP3, PP4, PP5, BA1, BS1, BS2, BS3, BS4, BP1, BP2, BP3, BP4, BP5, BP6, BP7)

# Autosomes ----
autosomes_pivoted_bi_edited <- autosomes_pivoted_bi %>%
  left_join(variants, by = "variant") %>%
  rowwise() %>%
  mutate(transcript = strsplit(Name, split = "\\(")[[1]][1]) %>%
  mutate(coding_protein = strsplit(Name, split = ":")[[1]][2]) %>%
  mutate(coding = strsplit(coding_protein, split = " ")[[1]][1]) %>%
  mutate(protein = strsplit(coding_protein, split = " ")[[1]][2]) %>%
  ungroup() %>%
  mutate(protein = str_remove(str_remove(protein, "\\("), "\\)")) %>%
  mutate(dbsnp = ifelse(!is.na(`RS#.(dbSNP)`),glue("rs{`RS#.(dbSNP)`}"), NA )) %>%
  rename("ref" = "ReferenceAlleleVCF", "alt" = "AlternateAlleleVCF") %>%
  select(hpos_surdez, `./.`, `0/0`,`0/1`, `1/1`, gene_names, transcript, coding, protein, variant, ref, alt, dbsnp, highestclinsig, hwep, qc_call_rate, var_qc_dp)

autosomes_pivoted_multi_edited <- autosomes_pivoted_multi %>%
  left_join({variants %>% select(variant, gene_names)}, by = "variant") %>%
  rowwise() %>%
  mutate(transcript = strsplit(Name, split = "\\(")[[1]][1]) %>%
  mutate(coding_protein = strsplit(Name, split = ":")[[1]][2]) %>%
  mutate(coding = strsplit(coding_protein, split = " ")[[1]][1]) %>%
  mutate(protein = strsplit(coding_protein, split = " ")[[1]][2]) %>%
  ungroup() %>%
  mutate(protein = str_remove(str_remove(protein, "\\("), "\\)")) %>%
  mutate(dbsnp = ifelse(!is.na(`RS#.(dbSNP)`),glue("rs{`RS#.(dbSNP)`}"), NA )) %>%
  # filter(!is.na(highestclinsig)) %>%
  select(hpos_surdez, `./._new`, `0/0_new`,`0/1`, `1/1`, gene_names, transcript, coding, protein, variant, ref, alt, dbsnp, highestclinsig, hwep, qc_call_rate, var_qc_dp) %>%
  distinct() %>%
  rename(`0/0` = `0/0_new`, `./.` = `./._new`)

autosomes_complete <- bind_rows(autosomes_pivoted_bi_edited, autosomes_pivoted_multi_edited) %>%
  pivot_wider(id_cols = c("variant", "highestclinsig", "ref", "alt", "dbsnp", "gene_names", "transcript", "coding", "protein", "var_qc_dp"),
              names_from = "hpos_surdez",
              values_from = c("./.", "0/0", "0/1", "1/1",  "hwep", "qc_call_rate")) %>%
  
  mutate(across(starts_with("qc_call_rate"), ~ replace_na(.x, "fail"))) %>%
  
  mutate(heterozygote_freq_yes = `0/1_yes` / (surdez_sample_size - `./._yes`)) %>%
  mutate(heterozygote_freq_no = `0/1_no` / (nao_surdez_sample_size - `./._no`)) %>%
  
  mutate(alt_freq_yes = ((2 * `1/1_yes`) + `0/1_yes`) / (2 * (surdez_sample_size - `./._yes`)) ) %>%
  mutate(alt_freq_no = ((2 * `1/1_no`) + `0/1_no`) / (2 * (nao_surdez_sample_size - `./._no`)) ) %>%
  
  mutate(missing_overall = `./._no` + `./._yes`) %>%
  mutate(het_overall = `0/1_no` +`0/1_yes`) %>%
  mutate(homalt_overall = `1/1_no` + `1/1_yes`) %>%

  mutate(hetfreq_overall = (het_overall) / (surdez_sample_size + nao_surdez_sample_size - missing_overall)) %>%

  mutate(alt_freq_overall = ((2 * homalt_overall) + het_overall) / (2 * (surdez_sample_size + nao_surdez_sample_size - missing_overall)) ) %>%
  
  mutate(variant2 = variant) %>%
  
  separate(variant2, into = c("Chr", "Start"), sep = ":", convert = TRUE) %>%

  left_join(annovar, by = c("Chr" = "Chr", "Start" = "Start", "ref" = "Ref", "alt" = "Alt")) %>%
  relocate(gene_names,
           variant,
           highestclinsig,
           CLNSIGCONF,
           # dbsnp,
           ref,
           alt,
           transcript,
           coding,
           protein,
           alt_freq_yes,
           heterozygote_freq_yes,
           `0/0_yes`,
           `0/1_yes`,
           `1/1_yes`, 
           `./._yes`, 
           qc_call_rate_yes,
           hwep_yes, 
           alt_freq_no,
           heterozygote_freq_no, 
           `0/0_no`, 
           `0/1_no`, 
           `1/1_no`, 
           `./._no`, 
           qc_call_rate_no, 
           hwep_no, 
           alt_freq_overall,
           hetfreq_overall, 
           `het_overall`, 
           `homalt_overall`, 
           `missing_overall`, 
           var_qc_dp) %>%
  mutate(tmp = str_remove(variant, "chr")) %>%
  separate(tmp, into = c("chr_n", "position"), convert = TRUE) %>%
  arrange(chr_n, position) %>%
  select(-c("Chr","Start", "chr_n", "position")) %>%
  rename(
    "Gene" = "gene_names",
    "Chromosome position" = "variant",
    "ACMG classification" = "highestclinsig",
    "CLNSIGCONF" = "CLNSIGCONF",
    # "dbSNP" = "dbsnp",
    "Reference allele" = "ref",
    "Alternate allele" = "alt",
    "Transcript" = "transcript",
    "cDNA" = "coding",
    "Protein" = "protein",
    
    "Alternate allele frequency (hearing loss group)" = "alt_freq_yes",
    "Heterozygote genotypic frequency (hearing loss group)" = "heterozygote_freq_yes",
    "N Reference homozygotes (0/0) (hearing loss group)" = "0/0_yes",
    "N Heterozygotes (0/1) (hearing loss group)" = "0/1_yes",
    "N Alternative homozygotes (1/1) (hearing loss group)" = "1/1_yes",
    "Not called (hearing loss group)" = "./._yes",
    "Call rate (hearing loss group)" = "qc_call_rate_yes",
    "HWE test p-value (hearing loss group)" = "hwep_yes",
    
    "Alternate allele frequency (no hearing loss group)" = "alt_freq_no",
    "Heterozygote genotypic frequency (no hearing loss group)" = "heterozygote_freq_no",
    "N Reference homozygotes (0/0) (no hearing loss group)" = "0/0_no",
    "N Heterozygotes (0/1) (no hearing loss group)" = "0/1_no",
    "N Alternative homozygotes (1/1) (no hearing loss group)" = "1/1_no",
    "Not called (no hearing loss group)" = "./._no",
    "Call rate (no hearing loss group)" = "qc_call_rate_no",
    "HWE test p-value (no hearing loss group)" = "hwep_no",
    
    "Alternate allele frequency (overall)" = "alt_freq_overall",
    "Heterozygote genotypic frequency (overall)" = "hetfreq_overall",
    "N Heterozygotes (0/1) (overall)" = "het_overall",
    "N Alternative homozygotes (1/1) (overall)" = "homalt_overall",
    "Not called (overall)" = "missing_overall",
    
    "gnomAD frequency (all)" = "gnomAD_genome_ALL",
    "ABRaOM frequency" = "ABRaOM_Frequencies",
    "REVEL score" = "REVEL_score",
    "REVEL rankscore" = "REVEL_rankscore",
    "SIFT score" = "SIFT_score",
    "SIFT pred" = "SIFT_pred",
    "SIFT converted rankscore" = "SIFT_converted_rankscore",
    "Polyphen2 HDIV score" = "Polyphen2_HDIV_score",
    "Polyphen2 HDIV rankscore" = "Polyphen2_HDIV_rankscore",
    "Polyphen2 HDIV pred" = "Polyphen2_HDIV_pred",
    "Polyphen2 HVAR score" = "Polyphen2_HVAR_score",
    "Polyphen2 HVAR rankscore" = "Polyphen2_HVAR_rankscore",
    "Polyphen2 HVAR pred" = "Polyphen2_HVAR_pred",
    
    "Variant QC" = "var_qc_dp"
  )
    
autosomes_filtered <- autosomes_complete %>%
  filter(`Variant QC` == "pass") %>%
  filter(`Call rate (hearing loss group)` == "pass" & `Call rate (no hearing loss group)` == "pass") %>%
  filter(`Alternate allele frequency (overall)` > 0) %>%
  filter(!is.na(Transcript)) %>%
  distinct()

# Chromosome X ----
chrx_pivoted_complete <- 
  chrx_pivoted %>%
  left_join(variants, by = "variant") %>%
  rowwise() %>%
  mutate(transcript = strsplit(Name, split = "\\(")[[1]][1]) %>%
  mutate(coding_protein = strsplit(Name, split = ":")[[1]][2]) %>%
  mutate(coding = strsplit(coding_protein, split = " ")[[1]][1]) %>%
  mutate(protein = strsplit(coding_protein, split = " ")[[1]][2]) %>%
  ungroup() %>%
  mutate(protein = str_remove(str_remove(protein, "\\("), "\\)")) %>%
  mutate(dbsnp = ifelse(!is.na(`RS#.(dbSNP)`),glue("rs{`RS#.(dbSNP)`}"), NA )) %>%
  rename("ref" = "ReferenceAlleleVCF", "alt" = "AlternateAlleleVCF") %>%
  select(hpos_surdez, `F_./.`, `F_0/0`,`F_0/1`, `F_1/1`, M_0, M_1, alt1_freq, heterozygote_freq, hemizygote_freq, gene_names, transcript, coding, protein, variant, ref, alt, dbsnp, highestclinsig, hwep, qc_call_rate, var_qc_dp) %>%
  pivot_wider(id_cols = c("variant", "highestclinsig", "ref", "alt", "dbsnp", "gene_names", "transcript", "coding", "protein", "var_qc_dp"),
              names_from = "hpos_surdez",
              values_from = c("F_./.", "F_0/0", "F_0/1", "F_1/1", "M_0", "M_1", "alt1_freq", "heterozygote_freq", "hemizygote_freq", "hwep", "qc_call_rate")) %>%
  
  mutate(M_._no = 0, M_._yes = 0) %>%
  mutate(missing_females_overall = `F_./._no` + `F_./._yes`) %>%
  mutate(missing_males_overall = M_._no + M_._yes) %>%
  mutate(homref_females_overall = `F_0/0_no` + `F_0/0_yes`) %>%
  mutate(het_females_overall = `F_0/1_no` +`F_0/1_yes`) %>%
  mutate(homalt_females_overall = `F_1/1_no` + `F_1/1_yes`) %>%
  mutate(hemiz_ref_males_overall = M_0_no + M_0_yes) %>%
  mutate(hemiz_alt_males_overall = M_1_no + M_1_yes) %>%
  mutate(hetfreq_females_overall = (het_females_overall) / (homref_females_overall + homalt_females_overall)) %>%
  mutate(hemiz_freq_overall = (hemiz_alt_males_overall) / (hemiz_alt_males_overall + hemiz_ref_males_overall)) %>%
  mutate(alt_freq_overall = ( 2 * (homalt_females_overall) + het_females_overall) / 2 * (homref_females_overall + het_females_overall + homalt_females_overall)) %>%
  
  mutate(variant2 = variant) %>%
  
  separate(variant2, into = c("Chr", "Start"), sep = ":", convert = TRUE) %>%
  left_join(annovar, by = c("dbsnp" = "dbsnp", "Chr" = "Chr", "Start" = "Start", "ref" = "Ref", "alt" = "Alt")) %>%
  
  relocate(gene_names, 
           variant, 
           highestclinsig, 
           CLNSIGCONF,
           dbsnp,
           ref, 
           alt, 
           transcript, 
           coding, 
           protein, 
           
           alt1_freq_yes,
           heterozygote_freq_yes, 
           hemizygote_freq_yes,
           `F_0/0_yes`,
           `F_0/1_yes`,
           `F_1/1_yes`,
           `F_./._yes`,
           M_0_yes,
           M_1_yes,
           M_._yes,
           qc_call_rate_yes,
           hwep_yes,
           
           alt1_freq_no,
           heterozygote_freq_no, 
           hemizygote_freq_no,
           `F_0/0_no`,
           `F_0/1_no`,
           `F_1/1_no`,
           `F_./._no`, 
           M_0_no,
           M_1_no,
           M_._no,
           qc_call_rate_no,
           hwep_no,
           
           alt_freq_overall,
           hetfreq_females_overall,
           hemiz_freq_overall,
           homref_females_overall,
           het_females_overall,
           homalt_females_overall,
           missing_females_overall,
           hemiz_ref_males_overall,
           hemiz_alt_males_overall,
           missing_males_overall,
           
           var_qc_dp
           
  ) %>%

    rename(
      "Gene" = "gene_names",
      "Chromosome position" = "variant",
      "ACMG classification" = "highestclinsig",
      "dbSNP" = "dbsnp",
      "Reference allele" = "ref",
      "Alternate allele" = "alt",
      "Transcript" = "transcript",
      "cDNA" = "coding",
      "Protein" = "protein",
      
      "Alternate allele frequency (hearing loss group)" = "alt1_freq_yes",
      "Females Heterozygote genotypic frequency (hearing loss group)" = "heterozygote_freq_yes",
      "Males Hemizygote frequency (hearing loss group)" = "hemizygote_freq_yes",
      "N Females reference homozygotes (0/0) (hearing loss group)" = "F_0/0_yes",
      "N Females heterozygotes (0/1) (hearing loss group)" = "F_0/1_yes",
      "N Females alternative homozygotes (1/1) (hearing loss group)" = "F_1/1_yes",
      "Not called females (hearing loss group)" = "F_./._yes",
      "N Males reference hemizygotes (0) (hearing loss group)" = "M_0_yes",
      "N Males alternate hemizygotes (1) (hearing loss group)" = "M_1_yes",
      "Not called males (hearing loss group)" = "M_._yes",
      "Call rate (hearing loss group)" = "qc_call_rate_yes",
      "HWE test p-value (hearing loss group)" = "hwep_yes",
      
      "Alternate allele frequency (no hearing loss group)" = "alt1_freq_no",
      "Females Heterozygote genotypic frequency (no hearing loss group)" = "heterozygote_freq_no",
      "Males Hemizygote frequency (no hearing loss group)" = "hemizygote_freq_no",
      "N Females reference homozygotes (0/0) (no hearing loss group)" = "F_0/0_no",
      "N Females heterozygotes (0/1) (no hearing loss group)" = "F_0/1_no",
      "N Females alternative homozygotes (1/1) (no hearing loss group)" = "F_1/1_no",
      "Not called females (no hearing loss group)" = "F_./._no",
      "N Males reference hemizygotes (0) (no hearing loss group)" = "M_0_no",
      "N Males alternate hemizygotes (1) (no hearing loss group)" = "M_1_no",
      "Not called males (no hearing loss group)" = "M_._no",
      "Call rate (no hearing loss group)" = "qc_call_rate_no",
      "HWE test p-value (no hearing loss group)" = "hwep_no",
      
      "Alternate allele frequency (overall)" = "alt_freq_overall",
      "Females Heterozygote genotypic frequency (overall)" = "hetfreq_females_overall",
      "Males Hemizygote frequency (overall)" = "hemiz_freq_overall",
      "N Females reference homozygotes (0/0) (overall)" = "homref_females_overall",
      "N Females heterozygotes (0/1) (overall)" = "het_females_overall",
      "N Females alternative homozygotes (1/1) (overall)" = "homalt_females_overall",
      "Not called females (overall)" = "homalt_females_overall",
      "N Males reference hemizygotes (0) (overall)" = "hemiz_ref_males_overall",
      "N Males alternate hemizygotes (1) (overall)" = "hemiz_alt_males_overall",
      "Not called males (overall)" = "missing_males_overall",
      
      "gnomAD frequency (all)" = "gnomAD_genome_ALL",
      "ABRaOM frequency" = "ABRaOM_Frequencies",
      "REVEL score" = "REVEL_score",
      "REVEL rankscore" = "REVEL_rankscore",
      "SIFT score" = "SIFT_score",
      "SIFT pred" = "SIFT_pred",
      "SIFT converted rankscore" = "SIFT_converted_rankscore",
      "Polyphen2 HDIV score" = "Polyphen2_HDIV_score",
      "Polyphen2 HDIV rankscore" = "Polyphen2_HDIV_rankscore",
      "Polyphen2 HDIV pred" = "Polyphen2_HDIV_pred",
      "Polyphen2 HVAR score" = "Polyphen2_HVAR_score",
      "Polyphen2 HVAR rankscore" = "Polyphen2_HVAR_rankscore",
      "Polyphen2 HVAR pred" = "Polyphen2_HVAR_pred",
      
      "Variant QC" = "var_qc_dp"
  ) %>%
  select(-c("Chr", "Start"))
  
chrx_pivoted_filtered <- chrx_pivoted_complete %>%
  filter(`Variant QC` == "pass") %>%
  filter(`Call rate (hearing loss group)` == "pass" & `Call rate (no hearing loss group)` == "pass") %>%
  filter(`Alternate allele frequency (overall)` > 0 | `Males Hemizygote frequency (overall)` > 0) %>%
  filter(!is.na(Transcript))
  
# CNV e SV: teste de associação genética para selecionar possíveis variantes novas de interesse

mt_info <- read.xlsx(here("data_tables", "gene_names.xlsx"), sheet = "mt") %>%
  select(variant, position, ref, alt)

chrm_pivoted_edit <- chrm_pivoted %>% 
  mutate(variant = str_replace(variant, "chrM", "chrMT")) %>%
  left_join(mt_info, by = "variant") %>%
  mutate(cDNA = glue("{position}{ref}>{alt}")) %>%
  pivot_wider(id_cols = c("gene_names", "variant", "ref", "alt", "cDNA", "var_qc_dp", "chr_n", "position"),
              names_from = "hpos_surdez",
              values_from = c("./.", "0/0", "call_rate")) %>%
  relocate(`0/0_yes`, `./._yes`, call_rate_yes, `0/0_no`, `./._no`, call_rate_no, .after = position) %>%
  relocate(var_qc_dp, .after = call_rate_no) %>%
  rename(
    "Chromosome position" = "variant",
    "Reference allele" = "ref",
    "Alternate allele" = "alt",
    
    "Individuals with reference calls (0/0) (hearing loss group)" = "0/0_yes",
    "Individuals with missing calls (hearing loss group)" = `./._yes`,
    "Call rate (hearing loss group)" = "call_rate_yes",
    
    "Individuals with reference calls (0/0) (no hearing loss group)" = "0/0_no",
    "Individuals with missing calls (no hearing loss group)" = `./._no`,
    "Call rate (no hearing loss group)" = "call_rate_no",
    
    "Variant QC" = "var_qc_dp"
  )
  
wb <- createWorkbook()

addWorksheet(wb, sheetName = "Autosomes (filt)", gridLines = TRUE)

writeData(wb, sheet = "Autosomes (filt)", x = autosomes_filtered)

# addWorksheet(wb, sheetName = "Autosomes (orig)", gridLines = TRUE)
# 
# writeData(wb, sheet = "Autosomes (orig)", x = autosomes_complete)

addWorksheet(wb, sheetName = "Chr X (filt)", gridLines = TRUE)

writeData(wb, sheet = "Chr X (filt)", x = chrx_pivoted_filtered)

# addWorksheet(wb, sheetName = "Chr X (orig)", gridLines = TRUE)
# 
# writeData(wb, sheet = "Chr X (orig)", x = chrx_pivoted_complete)

addWorksheet(wb, sheetName = "Mitochondria", gridLines = TRUE)

writeData(wb, sheet = "Mitochondria", x = chrm_pivoted_edit)

saveWorkbook(wb, here("output_tables","ClinVar variants.xlsx"), overwrite = TRUE)

