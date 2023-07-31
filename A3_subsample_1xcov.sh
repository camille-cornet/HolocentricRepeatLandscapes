#!/bin/bash

### We want 0.1X coverage
### Genome size estimates with flow cytometry (1 pg = 978 Mbp)
### Subsample using a different nb of read pairs for each species

#SBATCH -p normal.168h
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -J A3_subsample_1xcov
#SBATCH -o /home/kay/repeatexplorer/carex/A3_subsample_1xcov/A3_subsample_1xcov.out
#SBATCH --mail-user=camille.cornet@unine.ch
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

### Set directories
INDIR=/data3/camille/repeatexplorer/carex/
mkdir /data3/camille/repeatexplorer/carex/A3_subsample_1xcov/
OUTDIR=/data3/camille/repeatexplorer/carex/A3_subsample_1xcov/

### Create a samples file with nb of read pairs to use
paste -d '\t' $INDIR/A1_concatenate_lanes/samples.txt $OUTDIR/nb_reads.txt \
> $OUTDIR/samples_sizes.txt

### Define the function
subsample_1xcov() {
  SAMPLE=$1
  READS=$2
  seqtk sample -s 10 $INDIR/A2_fastp_polyG_trimming/${SAMPLE}_R1_polyGtrimmed.fastq $READS | \
  gzip -c > $OUTDIR/${SAMPLE}_R1_subsampled1xcov_run1.fastq.gz
  seqtk sample -s 10 $INDIR/A2_fastp_polyG_trimming/${SAMPLE}_R2_polyGtrimmed.fastq $READS | \
  gzip -c > $OUTDIR/${SAMPLE}_R2_subsampled1xcov_run1.fastq.gz
  ### Do the subsampling 2x per individual
  seqtk sample -s 20 $INDIR/A2_fastp_polyG_trimming/${SAMPLE}_R1_polyGtrimmed.fastq $READS | \
  gzip -c > $OUTDIR/${SAMPLE}_R1_subsampled1xcov_run2.fastq.gz
  seqtk sample -s 20 $INDIR/A2_fastp_polyG_trimming/${SAMPLE}_R2_polyGtrimmed.fastq $READS | \
  gzip -c > $OUTDIR/${SAMPLE}_R2_subsampled1xcov_run2.fastq.gz
}
### Subsample to 1x coverage for each file in parallel
export OUTDIR INDIR
export -f subsample_1xcov
parallel --colsep '\t' 'subsample_1xcov {1} {2}' :::: $OUTDIR/samples_sizes.txt
