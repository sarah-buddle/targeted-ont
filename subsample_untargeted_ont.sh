## Subsampling untargeted ONT data ##

# Calculate total basecalled length
awk  'NR%4==2 {print length}' ${results}/${run}/HAC_v0.9.0/DNA/basecall/${filename_dna}.fastq | 
awk '{SUM+=$1}END{print SUM}' \
> ${results}/${run}/HAC_v0.9.0/DNA/basecall/total_length.txt

awk  'NR%4==2 {print length}' ${results}/${run}/HAC_v0.9.0/RNA/basecall/${filename_rna}.fastq | 
awk '{SUM+=$1}END{print SUM}' \
> ${results}/${run}/HAC_v0.9.0/RNA/basecall/total_length.txt

# Subsample as if 16 samples had been sequenced per flow cell
# 45245565112 = ((84,416,431,871 + 60,369,376,488)/2) * 10 /16 to account for different sample numbers
python3 $subsample_py \
${results}/${run}/HAC_v0.9.0/DNA/basecall/${filename_dna}.fastq \
${results}/${run}/HAC_v0.9.0/DNA/basecall/${filename_dna}_sub.fastq \
45245565112

python3 $subsample_py \
${results}/${run}/HAC_v0.9.0/RNA/basecall/${filename_rna}.fastq \
${results}/${run}/HAC_v0.9.0/RNA/basecall/${filename_rna}_sub.fastq \
45245565112

# Check lengths of output
awk  'NR%4==2 {print length}' ${results}/${run}/HAC_v0.9.0/DNA/basecall/${filename_dna}_sub.fastq | 
awk '{SUM+=$1}END{print SUM}' \
> ${results}/${run}/HAC_v0.9.0/DNA/basecall/total_length_sub.txt

awk  'NR%4==2 {print length}' ${results}/${run}/HAC_v0.9.0/RNA/basecall/${filename_rna}_sub.fastq | 
awk '{SUM+=$1}END{print SUM}' \
> ${results}/${run}/HAC_v0.9.0/RNA/basecall/total_length_sub.txt

