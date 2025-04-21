# targeted-ont

Contains the scripts associated with the article: <b>Nanopore-based hybridization capture approaches for targeted viral metagenomics and whole-genome sequencing</b>.

## Targeted viral metagenomics
Raw pod5 files from ONT sequencing were basecalled and demultiplexed using dorado.sh.

Read and base counts were performed using counts_basecall.sh and read_counts.sh.

The raw data was input into ONT's EPI2ME pipeline (epi2me.sh) and into nf-core taxprofiler (taxprofiler_preprocess_nanopore.sh). Kraken2 and Bracken were also run through taxprofiler (taxprofiler_kraken.sh).

Taxprofiler preprocessed reads were aligned to the viral reference genomes and coverage and depths calculated using minimap.sh.

EPI2ME outputs were linked with timings data and coverage using read_timings.sh.

Alignment and coverage plots were produced with assigned_bases.R and coverage.R using the outputs of minimap.sh.

Read timings plots were produced with timings_byvirus_twist_ont.R using the outputs of read_timings.sh.

Outputs of taxonomic classifiers were filtered using metathresholds.R and plotted using taxonomic_classifiers.R.

Plots for comparison of different basecalling and demultiplexing approaches were produced with demux_comparison.R using the outputs of counts_basecall.sh.

## Untargeted metagenomics

Previously obtained ONT data was basecalled, subsampled and demultiplexed using dorado_untargeted.sh and subsample_untargeted_ont.sh, and was then analysed using the same methods as for targeted.

## Viral WGS
Data previously basecalled using Minknow was demultiplexed and trimmed using sureselect_preprocess.sh and then aligned to viral reference genomes using minimap.sh.







