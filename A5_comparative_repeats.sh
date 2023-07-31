#!/bin/bash

# The script must be launched from /data3/camille/repeatexplorer/carex/A5_comparative_repeats!

#SBATCH -p normal.1000h
#SBATCH -c 16
#SBATCH --mem=512G
#SBATCH -J A5_comparative_repeats
#SBATCH -o /home/kay/repeatexplorer/carex/A5_comparative_repeats/A5_comparative_repeats.out
#SBATCH --mail-user=camille.cornet@unine.ch
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

### Set directories
INDIR=/data3/camille/repeatexplorer/carex/A4_individual_repeats/
mkdir /data3/camille/repeatexplorer/carex/A5_comparative_repeats/
OUTDIR=/data3/camille/repeatexplorer/carex/A5_comparative_repeats/

### Rename each dataset to add sample name (species and individual)
for SAMPLE in $(cat /data3/camille/repeatexplorer/carex/A1_concatenate_lanes/samples.txt)
do
  ### Add prefix and concate
  seqtk rename $INDIR/${SAMPLE}_run1_merged.fasta ${SAMPLE}_ \
  > $OUTDIR/${SAMPLE}_run1_merged_renamed.fasta
done

### Concatenate all the samples fasta files
cat $OUTDIR/*run1_merged_renamed.fasta > $OUTDIR/all_samples_final.fasta

### Run RepeatExplorer for run1
singularity exec -e --bind $PWD:/data3/camille/repeatexplorer/carex/A5_comparative_repeats/ \
/home/kay/repeatexplorer/repeatexplorer_singularity/repex_tarean \
seqclust -p -c 16 -r 500000000 -tax VIRIDIPLANTAE3.0 --paired --prefix_length 8 \
-v /data3/camille/repeatexplorer/carex/A5_comparative_repeats/ \
/data3/camille/repeatexplorer/carex/A5_comparative_repeats/all_samples_final.fasta

### Vizualisation
#/home/camille/bin/revis-master/plot_comparative_clustering_summary.R \
#--cluster_table=/home/kay/repeatexplorer/A5_comparative_repeats/CLUSTER_TABLE.csv \
#--comparative_counts=/home/kay/repeatexplorer/A5_comparative_repeats/COMPARATIVE_ANALYSIS_COUNTS.csv \
#--number_of_colors=10 -n \
#--output=/home/kay/repeatexplorer/A5_comparative_repeats/COMPARATIVE_ANALYSIS_COUNTS_visualization_GS_normalized.pdf
### Cannot plot in same script because need package r-optparse and compatibility issues with singularity
