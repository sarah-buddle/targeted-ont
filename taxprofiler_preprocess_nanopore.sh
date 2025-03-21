## Run taxprofiler preprocessing ##

model=HAC_v0.9.0
#model=SUP_v0.9.0

for sample in ${samples[@]}; do

    echo $run $sample
    date

    mkdir ${results}/taxprofiler
    mkdir ${results}/taxprofiler/samplesheets

    echo -e "sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta" > ${results}/taxprofiler/samplesheets/samplesheet_${run}_test1.csv

    echo -e $sample","$run",OXFORD_NANOPORE,"$results"/"$model"/trim/"$sample"_"$run".fq.gz,," \
    >> ${results}/taxprofiler/samplesheets/samplesheet_${run}_test1.csv

    nextflow run nf-core/taxprofiler -r dev -profile conda -w ${results}/work --max_cpus 1 \
    -c ${software}/nf-core-configs/custom_resources.conf \
    --input ${results}/taxprofiler/samplesheets/samplesheet_${run}_test1.csv \
    --databases /SAN/breuerlab/BTRU-scratch/sarah/results/virus_methods/databases/database_kraken_refseq_nucleotide_v2.csv \
    --outdir ${results}/taxprofiler \
    --perform_longread_hostremoval \
    --hostremoval_reference ${human_genome} \
    --save_hostremoval_index \
    --save_hostremoval_unmapped 

done
