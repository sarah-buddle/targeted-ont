## Plot taxnomic classifier outputs ##

library(tidyverse)

run_levels <- c("twist_ont_251124", "nanopore_270923",)

run_labels <- c("ONT + Twist CVRP",  "Untargeted ONT")

tool_levels <- c("epi2me_bracken", "kraken2", "bracken", "czid")

tool_labels <- c("EPI2ME\n(Bracken)", "Kraken2", "Bracken", "CZ ID*")

raw_report <- read.csv("metathreholds.csv")

report_thresholds <- raw_report %>%
  mutate(count = 1) %>%
  filter(rank == "species" & !is.na(type) & !(type%in% c("Human", "Archaea")) & !(name_speciesorhigher %in% c("Lambdavirus lambda", "Emesvirus zinderi"))) %>%
  group_by(run, dilution, repeat., tool, type, result_category, basecall, demux, volume) %>%
  summarise(n_species = sum(count)) %>%
  ungroup() %>%
  pivot_wider(names_from = "result_category", values_from = "n_species", values_fill = 0) %>%
  mutate(threshold = "Thresholds")

report_nothresholds <- report_thresholds %>%
  mutate(false_positive = false_positive + true_negative,
         true_positive = true_positive + false_negative) %>%
  mutate(true_negative = 0, false_negative = 0, threshold = "No thresholds")

report <- rbind(report_thresholds, report_nothresholds) %>%
  select(-true_negative, -false_negative) %>%
  group_by(run, dilution, tool, type, threshold) %>%
  summarise(false_positive = mean(false_positive), true_positive = mean(true_positive)) %>%
  ungroup() %>%
  mutate(run = factor(run, run_levels, run_labels),
         tool = factor(tool, tool_levels, tool_labels),
         dilution = as.character(dilution))

viruses <- report %>%
  dplyr::filter(type == "Virus") %>%
  dplyr::mutate(sensitivity = true_positive / 6)

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

report %>%
  #filter(threshold == "Thresholds") %>%
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

ggsave(paste0(outdir, "false_positive_species_thres.png"), height = 4, width = 8)


false_positive_viruses <- raw_report %>%
  filter(type == "Virus" & result_category == "false_positive" & !(name_speciesorhigher %in% c("Lambdavirus lambda", "Emesvirus zinderi")))

fp_viruses_list <- false_positive_viruses %>%
  select(name, taxid) %>%
  unique() %>%
  arrange()

write.csv(fp_viruses_list, paste0(outdir, "/false_positive_viruses.csv"),  quote = FALSE, row.names = FALSE)

# Manual addition of virus types here

fp_virus_types <- read.csv(paste0(outdir, "/false_positive_viruses_types.csv")) %>%
  mutate(taxid = as.character(taxid))

fp_viruses <- false_positive_viruses %>%
  ungroup() %>%
  left_join(fp_virus_types) %>%
  mutate(count = 1) %>%
  group_by(dnarna_pair, tool, name_speciesorhigher, species_taxid, run, dilution, repeat., virus_type) %>%
  summarise(count = max(count)) %>%
  ungroup %>%
  group_by(dnarna_pair, tool, run, dilution, virus_type) %>%
  summarise(n_species = sum(count)) %>%
  ungroup() %>%
  group_by(tool, run, dilution, virus_type) %>%
  summarise(n_species = mean(n_species)) %>%
  ungroup() %>%
  mutate(run = factor(run, run_levels, run_labels),
         tool = factor(tool, tool_levels, tool_labels),
         dilution = as.character(dilution))

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




