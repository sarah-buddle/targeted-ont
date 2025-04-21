## Run filtering using metathresholds on taxonomic classifier outputs ##

devtools::load_all()

library(metathresholds)
library(tidyverse)

options(scipen=999)

samples <- read.csv("samples.csv")

full_report_thres <- makeFullReport(samplesheet_filepath, taxonomizr_sql, db_filepath, thresholds_filepath, collapse_species = TRUE,
                                    keep_species_only = FALSE, positive_species_filepath = positive_species_filepath, virus_only = FALSE)

combined_dnarna_report_thres <- combineDNARNA(full_report_thres)

full_report_withinfo <- combined_dnarna_report_thres %>%
  mutate(sample = str_extract(sample_id, "^(.*?)(?=_illumina|_nanopore|_twist|_twist_ont)"),
         run = str_extract(sample_id, "(?:illumina|nanopore|twist|twist_ont).*$")) %>%
  dplyr::mutate(run = ifelse(grepl("hours", sample), paste0(run, "_", str_extract(sample, "[^_]*$")), run)) %>%
  dplyr::mutate(sample = ifelse(grepl("hours", sample), str_extract(sample, "(.*?)_[^_]*"), sample)) %>%
  dplyr::left_join(samples)

write.csv(full_report_withinfo, "metathresholds.csv",
          quote = FALSE, row.names = FALSE)
