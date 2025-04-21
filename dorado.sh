## Basecall and demultiplex Twist ONT data ##

# HAC demultiplex
$dorado basecaller hac $data \
--emit-fastq \
-o ${results}/basecall \
--min-qscore 9 \
--no-trim 

# SUP demultiplex
$dorado basecaller sup $data \
--emit-fastq \
-o ${results}/basecall \
--min-qscore 10 \
--no-trim

model=HAC_v0.9.0
#model=SUP_v0.9.0

mkdir ${results}/${model}/demux ${results}/${model}/trim ${results}/${model}/demux_singleend

# Both ends with trimming
$dorado demux ${results}/${model}/basecall/${filename}.fastq \
--emit-fastq \
-o ${results}/${model}/demux \
--kit-name TWIST-96A-UDI \
--barcode-both-ends

# Both ends without trimming
$dorado demux ${results}/basecall/${filename}.fastq \
--emit-fastq \
-o ${results}/demux_notrim \
--kit-name TWIST-96A-UDI \
--barcode-both-ends \
--no-trim

# Custom trimming step
for filepath in ${results}/${model}/demux/*.fastq; do

    file=$(basename $filepath)

    $dorado trim $filepath \
    --emit-fastq \
    --primer-sequences $primers \
    -k SPM-CUSTOM \
    > ${results}/${model}/trim/${file}

done

# Single end
$dorado demux ${results}/${model}/basecall/${filename}.fastq \
--emit-fastq \
-o ${results}/${model}/demux_singleend \
--kit-name TWIST-96A-UDI

# Custom trimming
for filepath in ${results}/${model}/demux_singleend/*.fastq; do

    file=$(basename $filepath)

    sample=$(echo $file | sed 's/\.fastq//' | sed 's/.*_//')

    $dorado trim $filepath \
    --emit-fastq \
    --primer-sequences $primers \
    -k SPM-CUSTOM \
    > ${results}/${model}/trim/${sample}_singleend_all_${run}.fq

    gzip ${results}/${model}/trim/${sample}_singleend_all_${run}.fq

done




