## Create coverage plot ##

library(tidyverse)

twist_ont_1fc <- read.delim("coverage_twist_ont_1fc.txt", strip.white = TRUE,
                   header = FALSE,
                   col.names = c("sample", "run", "species", "rname",
                                 "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq")) %>%
  mutate(run = "twist_ont_251124")

ont_1fc <- read.delim("coverage_ont_1fc.txt", strip.white = TRUE,
                   col.names = c("sample", "run", "species", "rname",
                                 "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"),
                   header = FALSE) %>%
  mutate(run = "nanopore_270923")

samples <- read.csv("samples.csv")

run_levels <- c("twist_ont_251124", "nanopore_270923")

run_labels <- c( "ONT + Twist CVRP", "Untargeted ONT")

run_colours <-setNames(c("#619CFF", "#E6AB02", "#E6AB02", "#F8766D", "#00BA38"),
                       run_labels)


demux_levels <- c("HAC_v0.9.0_dualend", "HAC_v0.9.0_singleend", "SUP_v0.9.0_dualend", "SUP_v0.9.0_singleend", "Illumina_Illumina")

demux_labels <- c("HAC Dual-end", "HAC Single-end", "SUP Dual-end", "SUP Single-end", "Illumina")


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



format_labels <- function(l) {

  if (l == 0) {

    label <- 0

  } else {

    label <- 6 * 10^l

  }

  parse(text=label)

}

format_labels <- Vectorize(format_labels)

coverage <- rbind(twist_ont_1fc, ont_1fc) %>%
  left_join(samples) %>%
  filter(!is.na(dnarna)) %>%
  mutate(species = gsub("_", " ", species)) %>%
  mutate(n_flow_cells = ifelse(run == "twist_ont_18251124", "2", "1")) %>%
  mutate(dilution = as.numeric(dilution),
         species = factor(species, msa_2008_levels, msa_2008_labels),
         run = factor(run, run_levels, run_labels)) %>%
  filter((species %in% c("Human\nbetaherpesvirus 5", "Human\nmastadenovirus F") & dnarna == "DNA" &
            run %in% c("Untargeted Illumina", "Untargeted ONT")) |
           (species %in% c("Human\northopneumovirus", "Influenza B\nvirus", "Mammalian\northoreovirus", "Zika\nvirus")
            & dnarna == "RNA" & run %in% c("Untargeted Illumina", "Untargeted ONT")) |
           (run %in% c("Illumina + Twist CVRP", "ONT + Twist CVRP", "ONT + Twist CVRP (2 flow cells)")) &
           !(species %in% c("Lambda\nphage", "MS2\nphage"))) %>%
  dplyr::mutate(length = endpos) %>%
  select(-c(startpos, endpos, covbases)) %>%
  group_by(sample, run, species, dilution, basecall, demux, dnarna, n_flow_cells, volume) %>%
  summarise(coverage = weighted.mean(coverage, length)) %>%
  ungroup() %>%
  arrange(run, sample, species) %>%
  group_by(run, dilution, species, basecall, demux, dnarna, n_flow_cells, volume) %>%
  summarise(min_coverage = min(coverage), max_coverage = max(coverage), mean_coverage = mean(coverage)) %>%
  ungroup() %>%
  unite(basecall_demux, basecall, demux, sep = "_", remove = FALSE) %>%
  mutate(basecall_demux = factor(basecall_demux, demux_levels, demux_labels))

coverage %>%
  filter(dilution != 0) %>%
  filter(demux %in% c("dualend") & basecall %in% c("HAC_v0.9.0")) %>%
  mutate(dilution = ifelse(dilution == 0, 6, dilution)) %>%
  ggplot2::ggplot(mapping = aes(x = log10(dilution/6), y = mean_coverage, color = run, linetype = n_flow_cells)) +
  geom_smooth(lwd = 0.8, se = FALSE) +
  facet_wrap(vars(species)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top",
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  scale_x_continuous(breaks = c(1,2,3,4), labels = format_labels) +
  scale_color_manual(values = run_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 100)) +
  geom_errorbar(aes(ymin = min_coverage, max = max_coverage), width = 0.1) +
  ylab("Percentage genome coverage") +
  xlab("Genome copies per ml (gc/ml)")

ggsave(paste0(outdir, "coverage.png"), height = 4, width = 6.5)