## Comparing SUP vs HAC and single vs dual end demultiplexing ##

library(tidyverse)

# Import read and base counts data
counts <- read.csv("counts.csv", strip.white = TRUE)

# Levels and labels
sample_list <- c("barcode33", "barcode34", "barcode35", "barcode36", "barcode37", "barcode38", "barcode39", "barcode40")

demux_levels <- c("HAC_v0.9.0_trim_bothends", "HAC_v0.9.0_trim_singleend", "SUP_v0.9.0_trim_bothends", "SUP_v0.9.0_trim_singleend")

demux_labels <- c("HAC Dual-end", "HAC Single-end", "SUP Dual-end", "SUP Single-end")

# Convert filenames to basecalling and demuc data
find_bc_demux <- function(sample) {

  if (grepl("singleend_sup", sample)){
    return("SUP_v0.9.0_trim_singleend")
  } else if (grepl("singleend", sample)) {
    return("HAC_v0.9.0_trim_singleend")
  } else if (grepl("sup", sample)) {
    return("SUP_v0.9.0_trim_bothends")
  } else {
    return("HAC_v0.9.0_trim_bothends")
  }

}

find_bc_demux <- Vectorize(find_bc_demux)

# Format counts table
counts_table <- counts %>%
  select(-N50) %>%
  filter(n_reads != 0 & trim == "trim") %>%
  mutate(basecall_demux = find_bc_demux(sample)) %>%
  select(-basecall, -trim) %>%
  mutate(sample = ifelse(grepl(paste(sample_list, collapse = "|"), sample), str_extract(sample, "^[^_]+(?=_)"), "Other\nbarcodes")) %>%
  group_by(sample, run, basecall_demux) %>%
  summarise(n_reads = sum(n_reads), total_length = sum(total_length)) %>%
  ungroup() %>%
  mutate(basecall_demux = factor(basecall_demux, demux_levels, demux_labels))

# Write formatted counts table to csv
counts_table %>% 
  mutate(sample = gsub("\n", " ", sample)) %>% 
  write.csv(paste0(outdir, "demux_comparison.csv"), quote = FALSE, row.names = FALSE)

# Plot read counts with different basecalling/demux combinations
counts_table %>%
  filter(sample != "Other\nbarcodes") %>%
  ggplot(aes(x = sample, y = n_reads/10^6, fill = basecall_demux)) +
  geom_col(position = "dodge2") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = "Number of reads (millions)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_brewer(palette = "Paired")

ggsave(paste0(outdir, "/demux_comparison.png"), height = 3.5, width = 7)

# Plot number of other barcodes remaining
counts_table %>%
  filter(sample == "Other\nbarcodes") %>%
  ggplot(aes(x = sample, y = n_reads, fill = basecall_demux)) +
  geom_col(position = "dodge2") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = "Number of reads") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_brewer(palette = "Paired")

ggsave(paste0(outdir, "/demux_comparison_other.png"), height = 3.5, width = 3)
