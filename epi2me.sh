# Run EPI2ME wf-metagenomics

for sample in ${samples[@]}; do

        $nextflow run epi2me-labs/wf-metagenomics \
        -profile singularity \
        -w ${results}/epi2me/work \
        --fastq ${results}/HAC_v0.9.0/trim/${sample}_${time}_${run}.fq.gz \
        --classifier kraken2 \
        --exclude_host $human_genome \
        --database $kraken_db \
        --taxonomy $taxdump \
        --ref2taxid $ref2taxid \
        --include_kraken2_assignments \
        --out_dir ${results}/epi2me

done