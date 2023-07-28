#!/bin/bash

# The script must be launched from /data3/camille/repeatexplorer/carex/A4_individual_repeats!

#SBATCH -p normal.168h
#SBATCH -c 32
#SBATCH --mem=256G
#SBATCH -J A4_individual_repeats
#SBATCH -o /home/kay/repeatexplorer/carex/A4_individual_repeats/A4_individual_repeats.out
#SBATCH --mail-user=camille.cornet@unine.ch
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

### Set directories
INDIR=/data3/camille/repeatexplorer/carex/A3_subsample_1xcov/
mkdir /data3/camille/repeatexplorer/carex/A4_individual_repeats/
OUTDIR=/data3/camille/repeatexplorer/carex/A4_individual_repeats/

### Define the function
individual_repeat_analysis() {
  SAMPLE=$1
  ### Make interleaved fasta file for run1
  seqtk mergepe $INDIR/${SAMPLE}_R1_subsampled1xcov_run1.fastq.gz \
  $INDIR/${SAMPLE}_R2_subsampled1xcov_run1.fastq.gz | \
  seqtk seq -A > $OUTDIR/${SAMPLE}_run1_merged.fasta
  ### Make interleaved fasta file for run2
  seqtk mergepe $INDIR/${SAMPLE}_R1_subsampled1xcov_run2.fastq.gz \
  $INDIR/${SAMPLE}_R2_subsampled1xcov_run2.fastq.gz | \
  seqtk seq -A > $OUTDIR/${SAMPLE}_run2_merged.fasta
  ### Run RepeatExplorer for run1
  singularity exec -e --bind $PWD:/data3/camille/repeatexplorer/carex/A4_individual_repeats/ \
  /home/kay/repeatexplorer/repeatexplorer_singularity/repex_tarean \
  seqclust -p -c 16 -r 120000000 -tax VIRIDIPLANTAE3.0 -v /data3/camille/repeatexplorer/carex/A4_individual_repeats/${SAMPLE}_RE_run1 \
  /data3/camille/repeatexplorer/carex/A4_individual_repeats/${SAMPLE}_run1_merged.fasta
  ### Run RepeatExplorer for run2
  singularity exec -e --bind $PWD:/data3/camille/repeatexplorer/carex/A4_individual_repeats/ \
  /home/kay/repeatexplorer/repeatexplorer_singularity/repex_tarean \
  seqclust -p -c 16 -r 120000000 -tax VIRIDIPLANTAE3.0 -v /data3/camille/repeatexplorer/carex/A4_individual_repeats/${SAMPLE}_RE_run2 \
  /data3/camille/repeatexplorer/carex/A4_individual_repeats/${SAMPLE}_run2_merged.fasta
}
### Do the individual based repeat analysis for each individual in parallel
export OUTDIR INDIR TMPDIR
export -f individual_repeat_analysis
parallel --colsep '\t' 'individual_repeat_analysis {1}' --tmpdir $TMPDIR :::: /data3/camille/repeatexplorer/carex/A1_concatenate_lanes/samples.txt
