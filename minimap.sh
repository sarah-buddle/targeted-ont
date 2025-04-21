## Alignment to reference genomes for WGS data ###

mkdir -p ${results}/minimap/${run}/${sample}/${species}

# Align to viral genome using minimap - taprofiler filtered data for Twist CVRP, raw data for agilent sureselect 
${minimap} -ax map-ont ${genome_dir}/${species}_genome.fasta  \
${results}/trim/${sample}_${run}.fq.gz \
> ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}.sam

# Extract mapped reads
${samtools} view -bF 4 -h ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}.sam |
${samtools} sort \
> ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam

# Delete sam file to save space
if test -s ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam; then
    rm ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}.sam
fi

# Remove PCR duplicates
$samtools markdup -@ 4 -r -l 500000 ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam \
${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam \
-f ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup_stats.txt

# Index bam file
${samtools} index ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam 

# Calculate depth
${samtools} depth ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam \
> ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_depth.txt

# Calculate coverage
${samtools} coverage ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam \
> ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup_coverage.txt

# Calculate depth (no deduplication)
${samtools} depth ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam \
> ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_nodedup_depth.txt

# Calculate coverage (no deduplication)
${samtools} coverage ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam \
> ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_nodedup_coverage.txt

# Print coverage values to file
samplerun=$(echo -e $sample "\t" $run "\t" $species)

awk -v n="$samplerun" 'BEGIN{FS=OFS="\t"}{print n OFS $0}' \
${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup_coverage.txt |
tail -n +2 \
>> ${results}/minimap/coverage.txt

# Make fastq file of aligned reads
${samtools} fastq ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam  \
> ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.fastq

# Calculate number and length of aligned reads
minimap_reads=$(awk 'NR % 4 == 1' ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.fastq |
wc -l | awk '{print $1}')

minimap_bp=$(awk 'NR % 4 == 2' ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.fastq | 
tr -d '\n' | tr -d '[:space:]' | wc -m)

# Print to file
echo -e "minimap\t"$run"\t"$sample"\t"$species"\t"$minimap_reads"\t"$minimap_bp \
>> ${results}/minimap/all_minimap.txt

# As above for non-deduplicated reads
${samtools} fastq ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam  \
> ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_nodedup.fastq

minimap_reads_nodedup=$(awk 'NR % 4 == 1' ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_nodedup.fastq |
wc -l | awk '{print $1}')

minimap_bp_nodedup=$(awk 'NR % 4 == 2' ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_nodedup.fastq | 
tr -d '\n' | tr -d '[:space:]' | wc -m)

echo -e "minimap\t"$run"\t"$sample"\t"$species"\t"$minimap_reads_nodedup"\t"$minimap_bp_nodedup \
>> ${results}/minimap/all_minimap_nodedup.txt

# Create alignment plot using R script
$rscript $alignment_r \
${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_depth.txt \
${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_depth.png
