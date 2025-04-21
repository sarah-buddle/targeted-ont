## Plot EPI2ME output over time ##

library(tidyverse)
library(lubridate)

# Sample data
samples <- read.csv("samples.csv")

# Levels and labels
run_levels <- c("twist_ont_251124", "twist_ont_18251124")

run_labels <- c("1 flow cell", "2 flow cells")

repeat_levels <- c("A", "B")

repeat_labels <- c("Repeat A", "Repeat B")

virus_levels <- c("Human_herpesvirus_5",
                  "Human_mastadenovirus_F",
                  "Human_respiratory_syncytial_virus",
                  "Influenza_B_virus",
                  "Reovirus_3",
                  "Zika_virus")

virus_labels <- c("Human betaherpesvirus 5\ndsDNA\n229,354 nt\n1 segment",
                  "Human mastadenovirus F\ndsDNA\n34,392 nt\n1 segment",
                  "Human orthopneumovirus\n-RNA\n15,228 nt\n1 segment",
                  "Influenza B virus\n-RNA\n18,527 nt\n8 segments",
                  "Mammalian orthoreovirus\ndsRNA\n23,416 nt\n10 segments",
                  "Zika virus\n+RNA\n19,952 nt\n1 segment")

dilution_levels <- c("60000", "6000", "600", "60")


viral_reads_all <- data.frame()

for (run in c("twist_ont_18251124", "twist_ont_251124")) {

  # Read in coverage/timings data
  run_table <- read.delim(paste0(indir, run, "/timings/coverage_timings.txt"), sep = " ", header = FALSE,
                     col.names = c("read_id", "sample", "run", "species", "read_number", "read_length", "rname",
                                   "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq", "start_time")) %>%
    arrange(run, sample, species, start_time) %>%
    filter(!grepl("_singleend", sample)) %>%
    mutate(read_length = as.numeric(read_length)) %>%
    mutate(sample = paste0(sample, "_all"))

  # Format and calculate cumulative frequency
  viral_reads <- run_table %>%
    dplyr::mutate(across(where(is.character), str_trim)) %>%
    dplyr::filter(!grepl("unclassified", sample)) %>%
    dplyr::rename(chr_length = endpos) %>%
    dplyr::select(-c(read_number, startpos, numreads, covbases, meandepth, meanbaseq, meanmapq)) %>%
    dplyr::group_by(read_id, sample, run, species, read_length, start_time) %>%
    dplyr::summarise(coverage = weighted.mean(coverage, chr_length)) %>%
    dplyr::mutate(time_minutes = start_time / 60,
                  time_hours = start_time / 3600) %>%
    dplyr::mutate(species = factor(species, virus_levels, virus_labels)) %>%
    dplyr::arrange(start_time) %>%
    dplyr::group_by(sample, run, species) %>%
    dplyr::mutate(total_bases = cumsum(as.numeric(read_length))) %>%
    dplyr::ungroup()

  # Fill in missing values
  max_time <- floor(max(viral_reads$time_minutes, na.rm = TRUE))

  minute_marks <- seq(0, max_time, by = 1)

  empty <- viral_reads %>%
    distinct(sample, run, species) %>%
    mutate(total_bases = 0, coverage = 0, time_wholeminutes = 0, time_minutes = 0, time_hours = 0,
           read_length = NA)

  viral_reads_clean <- viral_reads %>%
    arrange(time_minutes) %>%
    select(-c(read_id, start_time)) %>%
    group_by(sample, run, species, time_wholeminutes = floor(time_minutes)) %>%
    slice(n()) %>%
    dplyr::filter(time_wholeminutes %in% minute_marks) %>%
    ungroup() %>%
    complete(sample, run, species, time_wholeminutes = minute_marks) %>%
    group_by(sample, run, species) %>%
    dplyr::mutate(time_wholeminutes = time_wholeminutes + 1) %>%
    rbind(empty) %>%
    arrange(time_wholeminutes) %>%
    group_by(sample, run, species) %>%
    mutate(total_bases = ifelse(is.na(total_bases) & lag(total_bases) == 0,  0, total_bases),
           coverage = ifelse(is.na(coverage) & lag(coverage) == 0,  0, coverage)) %>%
    fill(everything(), .direction = "down") %>%
    arrange(sample, run, species, species) %>%
    dplyr::left_join(samples) %>%
    mutate(run = factor(run, run_levels, run_labels),
           repeat. = factor(repeat., repeat_levels, repeat_labels))

  viral_reads_all <- rbind(viral_reads_all, viral_reads_clean)

}

# Plot

# Individual bases over time
for (dilution_value in c(60000, 6000, 600, 60)) {

  viral_reads_all %>%
    filter(dilution == dilution_value) %>%
    ggplot(aes(x = time_wholeminutes/60, y = total_bases/10^3, color = run, linetype = repeat.)) +
    geom_line() +
    facet_wrap(vars(species), ncol = 1, scales = "free_y", strip.position = "right") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          panel.grid.major.x = element_line(colour = "lightgrey", linewidth = 0.3),
          legend.title = element_blank(),
          strip.text.y = element_text(angle = 0)) +
    labs(x = "Time (hours)",
         y = "Total kbases") +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(0,72), breaks = c(0,12,24,36,48,60,72),
                       minor_breaks = c(0,12,24,36,48,60,72)) +
    scale_y_continuous(expand = expansion(mult = c(0, .2))) +
    ggtitle(paste0("Genome copies per ml: ", dilution_value))

  ggsave(paste0(outdir, "viruses_by_time_", dilution_value, "gcperml.png"),
         height = 6, width = 5.5)

}

# Individual bases over time - 5 hours
for (dilution_value in c(60000, 6000, 600, 60, 0)) {

  viral_reads_all %>%
    filter(dilution == dilution_value) %>%
    filter(time_wholeminutes <= 5*60) %>%
    ggplot(aes(x = time_wholeminutes/60, y = total_bases/10^3, color = run, linetype = repeat.)) +
    geom_line() +
    facet_wrap(vars(species), ncol = 1, scales = "free_y", strip.position = "right") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          panel.grid.major.x = element_line(colour = "lightgrey", linewidth = 0.3),
          panel.grid.minor.x = element_line(colour = "lightgrey", linewidth = 0.1),
          legend.title = element_blank(),
          strip.text.y = element_text(angle = 0)) +
    labs(x = "Time (hours)",
         y = "Total kbases") +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(0,1,2,3,4,5), minor_breaks = c(0,0.5,1,2,3,4,5)) +
    scale_y_continuous(expand = expansion(mult = c(0, .2))) +
    ggtitle(paste0("Genome copies per ml: ", dilution_value))

  ggsave(paste0(outdir, "viruses_by_time_5hours", dilution_value, "gcperml.png"),
         height = 6, width = 5.5)

}

# Combined coverage over time
  viral_reads_all %>%
    group_by(dilution, run, species, time_wholeminutes) %>%
    summarise(total_bases = mean(total_bases), coverage = mean(coverage)) %>%
    mutate(dilution = as.character(dilution)) %>%
    mutate(dilution = factor(dilution, dilution_levels)) %>%
    filter(dilution != "0") %>%
    ggplot(aes(x = time_wholeminutes/60, y = coverage, color = dilution, linetype = run)) +
    geom_line() +
    facet_wrap(vars(species), ncol = 1, strip.position = "right") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          panel.grid.major.x = element_line(colour = "lightgrey", linewidth = 0.3),
          strip.text.y = element_text(angle = 0)) +
    labs(x = "Time (hours)",
         y = "Genome coverage (%)",
         color = "Genome copies\nper ml",
         linetype = "Flow cells") +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(0,72), breaks = c(0,12,24,36,48,60,72)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0,100))

  ggsave(paste0(outdir, "coverage_by_time.png"),
         height = 8, width = 7)

  viral_reads_all %>%
    group_by(dilution, run, species, time_wholeminutes) %>%
    summarise(total_bases = mean(total_bases), coverage = mean(coverage)) %>%
    mutate(dilution = as.character(dilution)) %>%
    mutate(dilution = factor(dilution, dilution_levels)) %>%
    filter(dilution != "0") %>%
    filter(time_wholeminutes <= 5*60) %>%
    ggplot(aes(x = time_wholeminutes/60, y = coverage, color = dilution, linetype = run)) +
    geom_line() +
    facet_wrap(vars(species), ncol = 1, strip.position = "right") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          panel.grid.major.x = element_line(colour = "lightgrey", linewidth = 0.3),
          panel.grid.minor.x = element_line(colour = "lightgrey", linewidth = 0.1),
          strip.text.y = element_text(angle = 0)) +
    labs(x = "Time (hours)",
         y = "Genome coverage (%)",
         color = "Genome copies\nper ml",
         linetype = "Flow cells") +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(0,1,2,3,4,5), minor_breaks = c(0,0.5,1,2,3,4,5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0,100))

  ggsave(paste0(outdir, "coverage_by_time_5hours.png"),
         height = 8, width = 7)
