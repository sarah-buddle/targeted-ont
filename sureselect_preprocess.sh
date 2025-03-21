## Preprocess WGS data ##

kit=SQK-RPB114-24

for sample in ${samples[@]}; do

    zcat ${fastq_pass}/${sample}/*.fastq.gz |
    gzip \
    > ${data}/${run}/${sample}_${run}.fq.gz

    zcat ${data}/${run}/${sample}_${run}.fq.gz |
    $dorado demux \
    --emit-fastq \
    -o ${results}/demux_notrim \
    --kit-name $kit \
    --no-trim

    zcat ${data}/${run}/${sample}_${run}.fq.gz |
    $dorado demux \
    --emit-fastq \
    -o ${results}/demux \
    --kit-name $kit

    $dorado trim ${results}/demux/unknown_run_id_${kit}_${sample}.fastq \
    --emit-fastq \
    --primer-sequences $primers \
    -k SPM-CUSTOM \
    > ${results}/trim/${sample}_${run}.fq 

    gzip ${results}/trim/${sample}_${run}.fq

done



