## Plot taxnomic classifier outputs ##

library(tidyverse)

run_levels <- c("twist_ont_251124", "nanopore_270923",)

run_labels <- c("ONT + Twist CVRP",  "Untargeted ONT")

tool_levels <- c("epi2me_bracken", "kraken2", "bracken", "czid")

tool_labels <- c("EPI2ME\n(Bracken)", "Kraken2", "Bracken", "CZ ID*")

raw_report <- read.csv("metathreholds.csv")

# Import results with thresholds
report_thresholds <- raw_report %>%
  filter(run %in% c("twist_ont_251124", "nanopore_270923") & tool %in% tool_levels) %>%
  mutate(count = 1) %>%
  filter(rank == "species" & !is.na(type) & !(type%in% c("Human", "Archaea")) & !(name_speciesorhigher %in% c("Lambdavirus lambda", "Emesvirus zinderi"))) %>%
  group_by(run, dilution, repeat., tool, type, result_category, basecall, demux, volume) %>%
  summarise(n_species = sum(count)) %>%
  ungroup() %>%
  pivot_wider(names_from = "result_category", values_from = "n_species", values_fill = 0) %>%
  mutate(threshold = "Thresholds")

# Results if no thresholds had been applied
report_nothresholds <- report_thresholds %>%
  mutate(false_positive = false_positive + true_negative,
         true_positive = true_positive + false_negative) %>%
  mutate(true_negative = 0, false_negative = 0, threshold = "No thresholds")

# Create empty grid for correct calculation of avreages between replicates
samples <- report_thresholds %>%
  select(run, dilution, repeat., basecall, demux, volume) %>%
  distinct()

empty <- expand.grid(run = unique(report_thresholds$run),
                     dilution = unique(report_thresholds$dilution),
                     tool = unique(report_thresholds$tool),
                     type = unique(report_thresholds$type),
                     threshold = c("Thresholds", "No thresholds")) %>%
  full_join(samples) %>%
  mutate(false_positive = 0, true_positive = 0)

# Combine report
report <- rbind(report_thresholds, report_nothresholds) %>%
  select(-true_negative, -false_negative) %>%
  rbind(empty) %>%
  arrange(run, dilution, repeat., tool, type, threshold, desc(false_positive), desc(true_positive)) %>%
  distinct(run, dilution, repeat., tool, type, threshold, .keep_all = TRUE) %>%
  group_by(run, dilution, tool, type, threshold) %>%
  summarise(false_positive = mean(false_positive), true_positive = mean(true_positive)) %>%
  ungroup() %>%
  mutate(run = factor(run, run_levels, run_labels),
         tool = factor(tool, tool_levels, tool_labels),
         dilution = as.character(dilution))

# Viruses only and calculate sensitivity
viruses <- report %>%
  dplyr::filter(type == "Virus") %>%
  dplyr::mutate(sensitivity = true_positive / 6)

# Sensitivity plot
viruses %>%
  ggplot(aes(x = dilution, y = sensitivity, fill = threshold)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.6) +
  facet_grid(cols = vars(tool), rows = vars(run)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top",
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ylab("Sensitivity") +
  xlab("Genome copies per ml (gc/ml)")

ggsave(paste0(outdir, "sensitivity.png"), height = 4, width = 8)

# False positive viruses before and after thresholds
viruses %>%
  ggplot(aes(x = dilution, y = false_positive, fill = threshold)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.6) +
  facet_grid(cols = vars(tool), rows = vars(run)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top",
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ylab("Number of false positive viral species") +
  xlab("Genome copies per ml (gc/ml)")

ggsave(paste0(outdir, "false_positive_viruses.png"), height = 4, width = 8)

# All false positive species before thresholds
report %>%
  filter(threshold == "No thresholds") %>%
  ggplot(aes(x = dilution, y = false_positive, fill = type)) +
  geom_col(width = 0.6) +
  facet_grid(cols = vars(tool), rows = vars(run)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top",
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ylab("Number of false positive species") +
  xlab("Genome copies per ml (gc/ml)") +
  scale_fill_brewer(palette = "Dark2")

ggsave(paste0(outdir, "false_positive_species.png"), height = 4, width = 8)

# All false positive species after thresholds
report %>%
  filter(threshold == "Thresholds") %>%
  ggplot(aes(x = dilution, y = false_positive, fill = type)) +
  geom_col(width = 0.6) +
  facet_grid(cols = vars(tool), rows = vars(run)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top",
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ylab("Number of false positive species") +
  xlab("Genome copies per ml (gc/ml)") +
  scale_fill_brewer(palette = "Dark2")

ggsave(paste0(outdir, "false_positive_species_thres.png"), height = 4, width = 8)

# Select false positive viruses after thresholds
false_positive_viruses <- raw_report %>%
  filter(type == "Virus" & result_category == "false_positive" & !(name_speciesorhigher %in% c("Lambdavirus lambda", "Emesvirus zinderi")))

# Make list of names
fp_viruses_list <- false_positive_viruses %>%
  select(name, taxid) %>%
  unique() %>%
  arrange()

# Write to csv and manually add types
write.csv(fp_viruses_list, paste0(outdir, "/false_positive_viruses.csv"),  quote = FALSE, row.names = FALSE)

fp_virus_types <- read.csv(paste0(outdir, "/false_positive_viruses_types.csv"))

# Calculate number of viruses of each type
fp_viruses1 <- false_positive_viruses %>%
  ungroup() %>%
  left_join(fp_virus_types) %>%
  mutate(count = 1) %>%
  group_by(dnarna_pair, tool, run, dilution, virus_type) %>%
  summarise(n_species = sum(count)) %>%
  ungroup()

# Fill in empty values to correctly calculate averages
samples <- false_positive_viruses %>%
  select(dnarna_pair, run, dilution) %>%
  distinct()

empty <- expand.grid(tool = unique(fp_viruses1$tool),
                     virus_type = unique(fp_viruses1$virus_type),
                     dnarna_pair = unique(fp_viruses1$dnarna_pair)) %>%
  mutate(n_species = 0) %>%
  left_join(samples, by = "dnarna_pair")

# Average between replicates and format
fp_viruses2 <- fp_viruses1 %>%
  rbind(empty) %>%
  arrange(dnarna_pair, tool, run, dilution, virus_type, desc(n_species)) %>%
  distinct(dnarna_pair, tool, run, dilution, virus_type, .keep_all = TRUE) %>%
  group_by(tool, run, dilution, virus_type) %>%
  summarise(n_species = mean(n_species)) %>%
  ungroup() %>%
  mutate(run = factor(run, run_levels, run_labels),
         tool = factor(tool, tool_levels, tool_labels),
         dilution = as.character(dilution))

# Plot false positive virus types
fp_viruses %>%
  ggplot(aes(x = dilution, y = n_species, fill = virus_type)) +
  geom_col(width = 0.6) +
  facet_grid(cols = vars(tool), rows = vars(run)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top",
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = c(0,2,4,6,8,10)) +
  ylab("Number of false positive species") +
  xlab("Genome copies per ml (gc/ml)") +
  scale_fill_brewer(palette = "Dark2")

ggsave(paste0(outdir, "false_positive_virus_types.png"), height = 4, width = 8)




