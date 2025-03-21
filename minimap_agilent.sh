## Alignment to reference genomes for WGS data ###

mkdir -p ${results}/minimap_raw/${run}/${sample}/${species}

${minimap} -ax map-ont ${genome_dir}/${species}_genome.fasta  \
${results}/trim/${sample}_${run}.fq.gz \
> ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}.sam

${samtools} view -bF 4 -h ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}.sam |
${samtools} sort \
> ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam

if test -s ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam; then
    rm ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}.sam
fi

$samtools markdup -@ 4 -r -l 500000 ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam \
${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam \
-f ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup_stats.txt

${samtools} index ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam 

${samtools} depth ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam \
> ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_depth.txt

${samtools} coverage ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam \
> ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup_coverage.txt

${samtools} depth ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam \
> ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_nodedup_depth.txt

${samtools} coverage ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam \
> ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_nodedup_coverage.txt

samplerun=$(echo -e $sample "\t" $run "\t" $species)

awk -v n="$samplerun" 'BEGIN{FS=OFS="\t"}{print n OFS $0}' \
${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup_coverage.txt |
tail -n +2 \
>> ${results}/minimap_raw/coverage.txt

${samtools} fastq ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam  \
> ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.fastq

minimap_reads=$(awk 'NR % 4 == 1' ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.fastq |
wc -l | awk '{print $1}')

minimap_bp=$(awk 'NR % 4 == 2' ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.fastq | 
tr -d '\n' | tr -d '[:space:]' | wc -m)

echo -e "minimap\t"$run"\t"$sample"\t"$species"\t"$minimap_reads"\t"$minimap_bp \
>> ${results}/minimap_raw/all_minimap.txt

${samtools} fastq ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam  \
> ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_nodedup.fastq

minimap_reads_nodedup=$(awk 'NR % 4 == 1' ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_nodedup.fastq |
wc -l | awk '{print $1}')

minimap_bp_nodedup=$(awk 'NR % 4 == 2' ${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_nodedup.fastq | 
tr -d '\n' | tr -d '[:space:]' | wc -m)

echo -e "minimap\t"$run"\t"$sample"\t"$species"\t"$minimap_reads_nodedup"\t"$minimap_bp_nodedup \
>> ${results}/minimap_raw/all_minimap_nodedup.txt

$rscript $alignment_r \
${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_depth.txt \
${results}/minimap_raw/${run}/${sample}/${species}/${run}_${sample}_${species}_depth.png

