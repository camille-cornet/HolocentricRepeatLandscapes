#!/bin/bash

#SBATCH -p normal.168h
#SBATCH -c 10
#SBATCH --mem=16G
#SBATCH -J A2_fastp_polyG_trimming
#SBATCH -o /home/kay/repeatexplorer/carex/A2_fastp_polyG_trimming/A2_fastp_polyG_trimming.out
#SBATCH --mail-user=camille.cornet@unine.ch
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

### Set directories
INDIR=/data3/camille/repeatexplorer/carex/A1_concatenate_lanes/
mkdir /data3/camille/repeatexplorer/carex/A2_fastp_polyG_trimming/
OUTDIR=/data3/camille/repeatexplorer/carex/A2_fastp_polyG_trimming/

### Function to parallelize
trim_polyG() {
  SAMPLE=$1
  ### We might as well be stringent on the length and quality required,
  ### since we subsample afterwards anyway (might as well subsample only good reads)
  fastp --trim_poly_g -l 120 -q 30 \
  -i $INDIR/${SAMPLE}_R1_concat.fastq -o $OUTDIR/${SAMPLE}_R1_polyGtrimmed.fastq \
  -I $INDIR/${SAMPLE}_R2_concat.fastq -O $OUTDIR/${SAMPLE}_R2_polyGtrimmed.fastq
}
### Trim polyG tail for each file in parallel
export INDIR OUTDIR
export -f trim_polyG
parallel --colsep '\t' 'trim_polyG {1}' :::: $INDIR/samples.txt
