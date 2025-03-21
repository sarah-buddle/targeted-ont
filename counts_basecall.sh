#!/bin/bash
#$ -S /bin/bash
#$ -l tmem=10G,h_vmem=10G
#$ -l h_rt=24:00:00
#$ -N counts_basecall
#$ -o /SAN/breuerlab/BTRU-scratch/sarah/cluster
#$ -e /SAN/breuerlab/BTRU-scratch/sarah/cluster

run=twist_ont_251124
results=/SAN/breuerlab/BTRU-scratch/sarah/results/twist_ont/${run}

# model=SUP_v0.9.0
# filename=calls_2025-01-20_T14-10-00

model=HAC_v0.9.0
filename=calls_2025-01-16_T09-20-34

cat ${results}/${model}/basecall/${filename}.fastq |
awk  'NR%4==2 {print length}' | 
awk '{SUM+=$1}END{print SUM}' \
> ${results}/${model}/basecall/total_length.txt

cat ${results}/${model}/basecall/${filename}.fastq |
awk  'NR%4==1 {print $0}' | 
wc -l \
> ${results}/${model}/basecall/total_reads.txt

## Untargeted ##

run=nanopore_270923
results=/SAN/breuerlab/BTRU-scratch/sarah/results/twist_ont/${run}

model=HAC_v0.9.0

dnarna=DNA
filename=calls_2025-01-22_T18-11-59_sub

dnarna=RNA
filename=calls_2025-01-24_T13-38-51_sub

cat ${results}/${model}/${dnarna}/basecall/${filename}.fastq |
awk  'NR%4==1 {print $0}' | 
wc -l \
> ${results}/${model}/${dnarna}/basecall/total_reads_sub.txt