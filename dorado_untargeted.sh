### Basecall and demultiplex untargeted ONT data ###

# Basecall (no demultiplexing)
$dorado basecaller hac $data \
--emit-fastq \
-o ${results}/basecall/HAC_v0.9.0/${dnarna} \
--min-qscore 9 \
--no-trim

kit=SQK-RPB114-24

# Demultiplex (with trimming)
mkdir ${results}/HAC_v0.9.0/DNA/demux ${results}/HAC_v0.9.0/RNA/demux

$dorado demux ${results}/HAC_v0.9.0/${dnarna}/basecall/${filename}.fastq \
--emit-fastq \
-o ${results}/HAC_v0.9.0/${dnarna}/demux \
--kit-name ${kit} \
--barcode-both-ends