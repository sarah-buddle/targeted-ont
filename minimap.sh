#!/bin/bash
#$ -S /bin/bash
#$ -l tmem=1G,h_vmem=1G
#$ -l h_rt=24:00:00
#$ -N minimap
#$ -o /SAN/breuerlab/BTRU-scratch/sarah/cluster
#$ -e /SAN/breuerlab/BTRU-scratch/sarah/cluster

run=twist_ont_251124
data=/SAN/breuerlab/BTRU
results=/SAN/breuerlab/BTRU-scratch/sarah/results/twist_ont/${run}
genome_dir=/SAN/breuerlab/BTRU-scratch/sarah/data/mock_community/msa_2008/Genomes
minimap=/share/apps/genomics/minimap2/minimap2
samtools=/share/apps/genomics/samtools-1.18/bin/samtools

species_list=(Human_herpesvirus_5 Human_mastadenovirus_F Human_respiratory_syncytial_virus Influenza_B_virus Reovirus_3 Zika_virus)

#species_list=(MS2_phage Lambda_phage)

# mkdir ${results}/minimap
# echo -e "sample\trun\tspecies\tchr\tposition\tref_base\tdepth\tread_bases\tquality" > ${results}/minimap/mpileup.txt
# echo -e "sample\trun\tspecies\tconsensus_base" > ${results}/minimap/consensus.txt

#samples=( barcode33_sub barcode34_sub barcode35_sub barcode36_sub barcode37_sub barcode38_sub barcode39_sub barcode40_sub )

# samples=( barcode33_all barcode34_all barcode35_all barcode36_all barcode37_all barcode38_all barcode39_all barcode40_all )

# samples=(barcode01_sub barcode02_sub barcode03_sub barcode04_sub barcode05_sub barcode06_sub barcode07_sub barcode08_sub barcode09_sub barcode10_sub \
# barcode11_sub barcode12_sub barcode13_sub barcode14_sub barcode15_sub barcode16_sub barcode17_sub barcode18_sub barcode19_sub )

#samples=( barcode33_singleend_sup_all barcode34_singleend_sup_all barcode35_singleend_sup_all barcode36_singleend_sup_all barcode37_singleend_sup_all barcode38_singleend_sup_all barcode39_singleend_sup_all barcode40_singleend_sup_all)

samples=( barcode33_sup_all barcode34_sup_all barcode35_sup_all barcode36_sup_all barcode37_sup_all barcode38_sup_all barcode39_sup_all barcode40_sup_all \
barcode33_singleend_all barcode34_singleend_all barcode35_singleend_all barcode36_singleend_all barcode37_singleend_all barcode38_singleend_all barcode39_singleend_all barcode40_singleend_all)


#for sample in ${samples[@]}; do
for sample in barcode33_singleend_all; do

    for species in ${species_list[@]}; do

        # Minimap ####
        mkdir -p ${results}/minimap/${run}/${sample}/${species}

        $minimap -d ${genome_dir}/${species}.mmi ${genome_dir}/${species}.fasta 

        # ${minimap} -ax map-ont ${genome_dir}/${species}.fasta   \
        # ${data}/${run}/${sample}_${run}.fq.gz \
        # > ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}.sam

        ${minimap} -ax map-ont ${genome_dir}/${species}_genome.fasta  \
        ${results}/taxprofiler/samtools/fastq/${sample}_${run}.unmapped_other.fastq.gz \
        > ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}.sam

        ${samtools} view -bF 4 -h ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}.sam |
        ${samtools} sort \
        > ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam

        if test -s ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam; then
            rm ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}.sam
        fi

        $samtools markdup -@ 4 -r -l 500000 ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam \
         ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam \
        -f ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup_stats.txt

        ${samtools} index ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam 

        ${samtools} depth ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam \
        > ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_depth.txt

        ${samtools} coverage ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam \
        > ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup_coverage.txt

         ${samtools} depth ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam \
        > ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_nodedup_depth.txt

        ${samtools} coverage ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam \
        > ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_nodedup_coverage.txt

        samplerun=$(echo -e $sample "\t" $run "\t" $species)

        awk -v n="$samplerun" 'BEGIN{FS=OFS="\t"}{print n OFS $0}' \
        ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup_coverage.txt |
        tail -n +2 \
        >> ${results}/minimap/coverage.txt

        ${samtools} fastq ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam  \
        > ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.fastq

        minimap_reads=$(awk 'NR % 4 == 1' ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.fastq |
        wc -l | awk '{print $1}')

        minimap_bp=$(awk 'NR % 4 == 2' ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.fastq | 
        tr -d '\n' | tr -d '[:space:]' | wc -m)

        echo -e "minimap\t"$run"\t"$sample"\t"$species"\t"$minimap_reads"\t"$minimap_bp \
        >> ${results}/minimap/all_minimap.txt

         ${samtools} fastq ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped.bam  \
        > ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_nodedup.fastq

        minimap_reads_nodedup=$(awk 'NR % 4 == 1' ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_nodedup.fastq |
        wc -l | awk '{print $1}')

        minimap_bp_nodedup=$(awk 'NR % 4 == 2' ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_nodedup.fastq | 
        tr -d '\n' | tr -d '[:space:]' | wc -m)
        
        echo -e "minimap\t"$run"\t"$sample"\t"$species"\t"$minimap_reads_nodedup"\t"$minimap_bp_nodedup \
        >> ${results}/minimap/all_minimap_nodedup.txt

        $samtools mpileup ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam  \
        --fasta-ref ${genome_dir}/${species}_genome.fasta -a \
        >  ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mpileup.txt

        awk -v n="$samplerun" 'BEGIN{FS=OFS="\t"}{print n OFS $0}' \
         ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mpileup.txt \
        >> ${results}/minimap/mpileup.txt

        $samtools consensus -f fasta -a \
        ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_mapped_dedup.bam \
        > ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_consensus.fasta

        grep -v '>' ${results}/minimap/${run}/${sample}/${species}/${run}_${sample}_${species}_consensus.fasta |
        sed 's/./\0\n/g' |
        sed '/^$/d' |
        awk -v n="$samplerun" 'BEGIN{FS=OFS="\t"}{print n OFS $0}' - \
        >> ${results}/minimap/consensus.txt

        done

done