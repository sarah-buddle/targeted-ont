## Read statistics ##

kit=TWIST-96A-UDI
filetype=trim
model=HAC_v0.9.0
bc_demux=singleend
time=all

for sample in ${samples[@]}; do

    mv ${results}/${model}/${filetype}/*_${kit}_${sample}.fastq \
    ${results}/${model}/${filetype}/${sample}_all_${run}.fq

    zcat ${results}/${model}/${filetype}/${sample}_${time}_${run}.fq.gz |
    awk  'NR%4==2 {print length}' | 
    sort -nr \
    > ${results}/lengths/${sample}_${time}_${run}_lengths.txt

    N_READS=$(wc -l ${results}/lengths/${sample}_${time}_${run}_lengths.txt | awk '{print $1}')

    TOTAL_LENGTH=$(awk '{SUM+=$1}END{print SUM}' ${results}/lengths/${sample}_${time}_${run}_lengths.txt)

    MEDIAN=$(cat ${results}/lengths/${sample}_${time}_${run}_lengths.txt |
    awk '{a[NR]=$1} END{if (NR%2==1) {print a[(NR+1)/2]} else {print (a[NR/2]+a[NR/2+1])/2}}')

    Q3=$(cat ${results}/lengths/${sample}_${time}_${run}_lengths.txt |
    awk '{a[NR]=$1} END{print a[int(NR*0.25)]}')

    Q1=$(cat ${results}/lengths/${sample}_${time}_${run}_lengths.txt |
    awk '{a[NR]=$1} END{print a[int(NR*0.75)]}')

    CUMULATIVE_LENGTH=0

    ### Initialize a variable to store the n50 value
    N50=0

    ### Read through the sorted contig lengths one by one
    while read LENGTH; do
        ### Increment the cumulative length by the length of the current contig
        CUMULATIVE_LENGTH=$((CUMULATIVE_LENGTH + LENGTH))

        ### If the cumulative length is greater than or equal to half of the total length, then we have found the n50 value
        if [ $CUMULATIVE_LENGTH -ge $((TOTAL_LENGTH / 2)) ]; then
        N50=$LENGTH
        break
        fi
    done < ${results}/lengths/${sample}_${time}_${run}_lengths.txt

    echo -e $sample"_"$time","$run","$model","$filetype","$N_READS","$TOTAL_LENGTH","$N50","$Q1","$MEDIAN","$Q3 \
    >> ${results}/counts.csv

done