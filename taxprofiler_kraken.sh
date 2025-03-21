## Run kraken2 through taxprofiler and bracken ##

for sample in ${samples[@]}; do

    echo $run $sample
    date

    echo -e "sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta" > ${results}/samplesheets/samplesheet_filtered_${run}_test1.csv

    echo -e $sample","$run",OXFORD_NANOPORE,"$results"/samtools/fastq/"$sample"_"$run".unmapped_other.fastq.gz,,," \
    >> ${results}/samplesheets/samplesheet_filtered_${run}_test1.csv

    nextflow run nf-core/taxprofiler -r dev -profile singularity -with-tower -w ${results}/work --max_cpus 1 \
    --input ${results}/samplesheets/samplesheet_filtered_${run}_test1.csv \
    --databases /SAN/breuerlab/BTRU-scratch/sarah/results/virus_methods/databases/database_kraken_refseq_nucleotide_v2.csv \
    --outdir ${results} \
    --run_kraken2 --kraken2_save_readclassifications \
    --run_bracken

    date

done

# Run bracken separately since it's not supported for ONT in taxprofiler
mamba activate kraken

# Threshold of 10 for Twist ONT and default 0 for untargeted ONT
for sample in ${samples[@]}; do

    bracken -d /SAN/breuerlab/BTRU-scratch2/sarah/software/kraken_refseq_db/refseq-2023-06-08-nucleotide-v2 \
    -i ${results}/kraken2/refseq-2023-06-08-nucleotide-v2/${sample}_se_${run}_refseq-2023-06-08-nucleotide-v2.kraken2.kraken2.report.txt \
    -o ${results}/bracken/refseq-2023-06-08-nucleotide-v2/${sample}_se_${run}_refseq-2023-06-08-nucleotide-v2.bracken.tsv \
    -t 10

done