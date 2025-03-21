library(tidyverse)

##### Levels and labels ####

run_levels <- c("twist_ont_251124", "nanopore_270923")

run_labels <- c( "ONT + Twist CVRP", "Untargeted ONT")

run_colours <-setNames(c("#619CFF", "#E6AB02", "#E6AB02", "#F8766D", "#00BA38"),
                       run_labels)

demux_levels <- c("HAC_v0.9.0_dualend", "HAC_v0.9.0_singleend", "SUP_v0.9.0_dualend", "SUP_v0.9.0_singleend")

demux_labels <- c("HAC Dual-end", "HAC Single-end", "SUP Dual-end", "SUP Single-end")


format_labels <- function(l) {

  if (l == 0) {

    label <- 0

  } else {

    label <- 6 * 10^l

  }

  parse(text=label)

}

format_labels <- Vectorize(format_labels)


msa_2008_levels <- c("Human betaherpesvirus 5", "Human herpesvirus 5", "Cytomegalovirus humanbeta5",
                     "Human mastadenovirus F",
                     "Human orthopneumovirus", "Orthopneumovirus hominis", "Human respiratory syncytial virus",
                     "Influenza B virus", "Betainfluenzavirus influenzae",
                     "Mammalian orthoreovirus", "Reovirus 3",
                     "Zika virus", "Orthoflavivirus zikaense",
                     "Lambda phage", "Escherichia virus Lambda", "Lambdavirus lambda",
                     "MS2 phage", "Escherichia virus MS2", "Emesvirus zinderi")

msa_2008_labels_oneline <- c("Human betaherpesvirus 5", "Human betaherpesvirus 5", "Human betaherpesvirus 5",
                             "Human mastadenovirus F",
                             "Human orthopneumovirus", "Human orthopneumovirus", "Human orthopneumovirus",
                             "Influenza B virus", "Influenza B virus",
                             "Mammalian orthoreovirus", "Mammalian orthoreovirus",
                             "Zika virus", "Zika virus",
                             "Lambda phage", "Lambda phage", "Lambda phage",
                             "MS2 phage", "MS2 phage", "MS2 phage")

msa_2008_labels <- c("Human\nbetaherpesvirus 5", "Human\nbetaherpesvirus 5", "Human\nbetaherpesvirus 5",
                     "Human\nmastadenovirus F",
                     "Human\northopneumovirus", "Human\northopneumovirus", "Human\northopneumovirus",
                     "Influenza B\nvirus", "Influenza B\nvirus",
                     "Mammalian\northoreovirus", "Mammalian\northoreovirus",
                     "Zika\nvirus", "Zika\nvirus",
                     "Lambda\nphage", "Lambda\nphage", "Lambda\nphage",
                     "MS2\nphage", "MS2\nphage", "MS2\nphage")


#### File inputs ####
samples <- read.csv("samples.csv")

#### BPM plot ####

twist_ont_1fc <- read.delim("all_minimap_twist_ont_1fc.txt",
                     sep = "\t", row.names = NULL,
                   header = FALSE, col.names = c("tool", "run", "sample", "species", "reads", "bases")) %>% 
  filter(!grepl("_sub", sample))

ont_1fc <- read.delim("all_minimap_ont_1fc.txt",
                   sep = "\t", row.names = NULL,
                   header = FALSE, col.names = c("tool", "run", "sample", "species", "reads", "bases"))

alignments <- rbind(twist_ont_1fc, ont_1fc) %>%
  mutate(bases = as.numeric(bases)) %>% 
  left_join(samples) %>%
  filter(!is.na(dnarna)) %>%
  mutate(species = gsub("_", " ", species)) %>%
  mutate(n_flow_cells = ifelse(run == "twist_ont_18251124", "2", "1")) %>%
  unite(basecall_demux, basecall, demux, sep = "_", remove = FALSE) %>%
  mutate(dilution = as.numeric(dilution),
         species = factor(species, msa_2008_levels, msa_2008_labels),
         run = factor(run, run_levels, run_labels),
         basecall_demux = factor(basecall_demux, demux_levels, demux_labels))) %>%
  filter((species %in% c("Human\nbetaherpesvirus 5", "Human\nmastadenovirus F") & dnarna == "DNA" &
            run %in% c("Untargeted Illumina", "Untargeted ONT")) |
           (species %in% c("Human\northopneumovirus", "Influenza B\nvirus", "Mammalian\northoreovirus", "Zika\nvirus")
            & dnarna == "RNA" & run %in% c("Untargeted Illumina", "Untargeted ONT")) |
           (run %in% c("Illumina + Twist CVRP", "ONT + Twist CVRP", "ONT + Twist CVRP (2 flow cells)")) &
           !(species %in% c("Lambda\nphage", "MS2\nphage"))) %>%
  group_by(run, species, dilution, dnarna, basecall_demux, n_flow_cells, volume) %>%
  summarise(min_bases = min(bases), max_bases = max(bases), mean_bases = mean(bases)) %>%
  ungroup() %>% 
  mutate(log10_bases = ifelse(mean_bases == 0, 0, log10(mean_bases)),
         log10_max_bases = ifelse(max_bases == 0, 0, log10(max_bases)),
         log10_min_bases = ifelse(min_bases == 0, 0, log10(min_bases))) %>%
  mutate(dilution = as.numeric(dilution)) %>%
  ungroup()

# Standard plot
alignments %>%
  filter(dilution != 0) %>%
  filter(basecall_demux %in% c("HAC Dual-end")) %>%
  mutate(dilution = ifelse(dilution == 0, 6, dilution)) %>%
  ggplot(aes(x = log10(dilution/6), y = log10_bases, color = run, linetype = n_flow_cells)) +
  facet_wrap(vars(species)) +
  geom_smooth(linewidth = 0.8, se = FALSE) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "top") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  scale_x_continuous(breaks = c(1,2,3,4), labels = format_labels) +
  scale_color_manual(values = run_colours) +
  xlab("Genome copies per ml (gc/ml)") +
  ylab("log10(Bases aligned to genome)") +
  geom_errorbar(aes(ymin = log10_min_bases, max = log10_max_bases), width = 0.1) +
  guides(fill = guide_legend(nrow = 2))

ggsave(paste0(outdir, "/bases.png"), height = 4, width = 6.5)

# Demux comparison table
basecall_demux <- alignments %>%
  filter(run == "ONT + Twist CVRP") %>%
  select(-c(min_bases, max_bases, log10_bases, log10_min_bases, log10_max_bases, dnarna)) %>%
  mutate(basecall_demux= gsub(" ", "_", basecall_demux)) %>%
  mutate(basecall_demux= gsub("-", "_", basecall_demux)) %>%
  pivot_wider(names_from = "basecall_demux", values_from = "mean_bases") %>%
  mutate(SUP_increase = (SUP_Dual_end - HAC_Dual_end) / HAC_Dual_end,
         SE_increase = (HAC_Single_end - HAC_Dual_end) / HAC_Dual_end,
         SUP_SE_increase = (SUP_Single_end - HAC_Dual_end) / HAC_Dual_end) %>%
  mutate(species = factor(species, msa_2008_labels, msa_2008_labels_oneline)) %>%
  arrange(species, desc(dilution)) %>%
  mutate_if(is.numeric, round, 2)

write.csv(basecall_demux, paste0(outdir, "/basecall_comparison_table.csv"),
          row.names = FALSE, quote = FALSE)

# PCR duplicates
twist_1fc_nodedup <- read.delim("/all_minimap_nodedup_twist_ont_1fc.txt",
                   sep = "\t", row.names = NULL,
                   header = FALSE, col.names = c("tool", "run", "sample", "species", "reads", "bases")) %>%
  mutate(tool = "minimap_nodedup")

ont_1fc_nodedup <- read.delim("all_minimap_nodedup_ont_1fc.txt",
                           sep = "\t", row.names = NULL,
                           header = FALSE, col.names = c("tool", "run", "sample", "species", "reads", "bases")) %>%
  mutate(tool = "minimap_nodedup")

duplicates <- rbind(twist_1fc, twist_1fc_nodedup, ont_1fc, ont_1fc_nodedup) %>%
  left_join(samples) %>%
  filter(!is.na(dnarna)) %>%
  mutate(species = gsub("_", " ", species)) %>%
  unite(basecall_demux, basecall, demux, sep = "_", remove = FALSE) %>%
  mutate(dilution = as.numeric(dilution),
         species = factor(species, msa_2008_levels, msa_2008_labels),
         run = factor(run, run_levels, run_labels),
         basecall_demux = factor(basecall_demux, demux_levels, demux_labels)) %>%
  filter((species %in% c("Human\nbetaherpesvirus 5", "Human\nmastadenovirus F") & dnarna == "DNA" &
            run %in% c("Untargeted Illumina", "Untargeted ONT")) |
           (species %in% c("Human\northopneumovirus", "Influenza B\nvirus", "Mammalian\northoreovirus", "Zika\nvirus")
            & dnarna == "RNA" & run %in% c("Untargeted Illumina", "Untargeted ONT")) |
           (run %in% c("Illumina + Twist CVRP", "ONT + Twist CVRP", "ONT + Twist CVRP (2 flow cells)")) &
           !(species %in% c("Lambda\nphage", "MS2\nphage"))) %>%
  group_by(tool, run, species, dilution, dnarna, basecall_demux) %>%
  summarise(mean_bases = mean(bases)) %>%
  ungroup() %>%
  pivot_wider(names_from = "tool", values_from = "mean_bases") %>%
  mutate(percent_duplicates = (minimap_nodedup - minimap)* 100 / minimap_nodedup) %>%
  mutate(percent_duplicates = ifelse(minimap_nodedup == 0, 0, percent_duplicates))

duplicates %>%
  ggplot(aes(x = log10(dilution/6), y = percent_duplicates, color = run)) +
  facet_wrap(vars(species)) +
  geom_smooth(linewidth = 0.8, se = FALSE) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "top") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  scale_x_continuous(breaks = c(1,2,3,4), labels = format_labels) +
  scale_color_manual(values = run_colours) +
  xlab("Genome copies per ml (gc/ml)") +
  ylab("Percentage duplicates")

ggsave(paste0(outdir, "/duplicates.png"), height = 4, width = 6.5)

