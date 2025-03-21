## Extract read timings data and coverage from EPI2ME results ##

declare -A taxons=( ["Human_mastadenovirus_F"]="130309" ["Human_herpesvirus_5"]="10359|3050295" \
["Human_respiratory_syncytial_virus"]="11250|3049954" ["Influenza_B_virus"]="518987|11520|2955465" \
["Zika_virus"]="64320|3048459" ["Reovirus_3"]="10886|538123|351073")

mkdir ${results}/timings

# Make compact sequencing summary file containing just read ids and start times
cut -f 5,10 ${seq_sum} |
tail -n +2 \
> ${results}/timings/seqsum_nohead_v1.txt

cut -f 5,10 ${seq_sum2} |
tail -n +2 \
>> ${results}/timings/seqsum_nohead_v1.txt

sort -k 1 ${results}/timings/seqsum_nohead_v1.txt \
> ${results}/timings/seqsum_nohead.txt

rm ${results}/timings/seqsum_nohead_v1.txt

echo -e "sample\trun\ttaxid\tread_id" > ${results}/timings/taxon_read_ids.txt

# Extract read IDs by taxid from EPI2ME output
for sample in ${samples[@]}; do

    for i in "${!taxons[@]}"; do 

        species=$i
        taxid=${taxons[$i]}

        awk -F'\t' -v s="$taxid" '$3 ~ s' ${results}/epi2me/kraken_reads_assignments/${sample}_all_${run}_lineages.kraken2.assignments.tsv |
        cut -f2 | sort -u \
        > ${results}/timings/${sample}_all_${run}_${species}_read_ids_v1.txt 

        samplerun=$(echo -e $sample "\t" $run "\t" $species)

        awk -v n="$samplerun"  'BEGIN{FS=OFS="\t"}{print n OFS $0}' ${results}/timings/${sample}_all_${run}_${species}_read_ids_v1.txt   \
        >> ${results}/timings/taxon_read_ids.txt

        rm ${results}/timings/${sample}_all_${run}_${species}_read_ids_v1.txt 

    done

done

# Join timings data to classified reads
tail -n +2 ${results}/timings/taxon_read_ids.txt |
grep "$run" |
sort -k 4 |
join - ${results}/timings/seqsum_nohead.txt -a 1 -1 4 -2 1 \
> ${results}/timings/taxon_read_ids_timings.txt

# Sort by timings
sort -k 5 -n ${results}/timings/taxon_read_ids_timings.txt \
> ${results}/timings/taxon_read_ids_timings_sorted.txt

# Add headings
echo -e "rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq" \
> ${results}/timings/coverage.txt

for sample in ${samples[@]}; do

     for i in "${!taxons[@]}"; do 

       species=$i
       taxid=${taxons[$i]}

        # Extract read ids belonging to that taxid
        awk -v sample="$sample" -v species="$species" '$2 == sample && $4 == species' \
        ${results}/timings/taxon_read_ids_timings_sorted.txt |
        sort -k 5 -n | awk '{print $1}' \
        > ${results}/timings/${sample}_all_${run}_${species}_read_ids.txt

        # Extract those reads from the fastq file
        zcat ${results}/${model}/trim/${sample}_all_${run}.fq.gz |
        awk 'NR % 4 == 1 {sub(/\t.*/, "");} {print}' |
        $seqtk subseq - ${results}/timings/${sample}_all_${run}_${species}_read_ids.txt \
        > ${results}/timings/${sample}_all_${run}_${species}_unsorted.fq

        # # Loop through each read ID and extract the corresponding read from the FASTQ file- so they are ordered by time sequenced
        
        while read -r id; do

            cat ${results}/timings/${sample}_all_${run}_${species}_unsorted.fq |
            grep "$id" -A 3 --no-group-separator \
            >> ${results}/timings/${sample}_all_${run}_${species}.fq
            
        done < ${results}/timings/${sample}_all_${run}_${species}_read_ids.txt

        # Count total reads assigned
        max_n_reads=$(wc -l ${results}/timings/${sample}_all_${run}_${species}_read_ids.txt | awk '{print $1}')

        # Loop through adding one additional read each time and calculating coverage
        for n_reads in $(seq 1 $max_n_reads); do

            # Calculate number of lines
            let n_lines=$n_reads*4

            # Extract reads, align to genome, calculate coverage with samtools
            head -${n_lines} ${results}/timings/${sample}_all_${run}_${species}.fq |
            $minimap -ax map-ont ${genome_dir}/${species}_genome.fasta - |
            $samtools sort |
            $samtools coverage - |
            tail -n +2 \
            > ${results}/timings/${sample}_all_${run}_${species}_coverage_temp.txt

            # Find corresponding read ID so we can link back to timings later
            read_id=$(head -${n_lines} ${results}/timings/${sample}_all_${run}_${species}.fq | tail -4 | head -1 | sed 's/@//')

            # Find corresponding read length after trimming
            read_length=$(head -${n_lines} ${results}/timings/${sample}_all_${run}_${species}.fq | tail -3 | head -1 | wc -c)

            # Combine into single file
            samplerun=$(echo -e $sample "\t" $run "\t" $species "\t" $n_reads "\t" $read_id "\t" $read_length)

            awk -v n="$samplerun" 'BEGIN{FS=OFS="\t"}{print n OFS $0}' \
            ${results}/timings/${sample}_all_${run}_${species}_coverage_temp.txt \
            >> ${results}/timings/coverage.txt

        done

    # Tidy files
    gzip ${results}/timings/${sample}_all_${run}_${species}.fq

    done

done

# Link output with timings data again
tail -n +2 ${results}/timings/coverage.txt |
sort -k 5 |
join - ${results}/timings/seqsum_nohead.txt -a 1 -1 5 -2 1 \
> ${results}/timings/coverage_timings.txt