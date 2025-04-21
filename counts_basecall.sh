#### Calculate length and bases for total basecalled file ####

cat ${results}/${model}/basecall/${filename}.fastq |
awk  'NR%4==2 {print length}' | 
awk '{SUM+=$1}END{print SUM}' \
> ${results}/${model}/basecall/total_length.txt

cat ${results}/${model}/basecall/${filename}.fastq |
awk  'NR%4==1 {print $0}' | 
wc -l \
> ${results}/${model}/basecall/total_reads.txt
