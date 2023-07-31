#!/bin/bash

#SBATCH -p normal.168h
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -J A1_concatenate_lanes
#SBATCH -o /home/kay/repeatexplorer/carex/A1_concatenate_lanes/A1_concatenate_lanes.out
#SBATCH --mail-user=camille.cornet@unine.ch
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

### Set directories
INDIR=/data3/kay
mkdir /data3/camille/repeatexplorer/carex/A1_concatenate_lanes/
OUTDIR=/data3/camille/repeatexplorer/carex/A1_concatenate_lanes/

### Set up files with all samples
find $INDIR -name "*.fastq.gz" | sort | rev | cut -d'_' -f9 | rev | uniq \
> $OUTDIR/inds.txt
paste -d '_' $OUTDIR/species.txt $OUTDIR/inds.txt > $OUTDIR/samples.txt

### Function to parallelize
concat() {
  SPE=$1
  IND=$2
  # Concatenate the fastq files for each lines for all individuals and species
  # (keeping forward and reverse in separate files)
  zcat $INDIR/BSSE*${IND}_*R1*.fastq.gz >> $OUTDIR/${SPE}_${IND}_R1_concat.fastq
  zcat $INDIR/BSSE*${IND}_*R2*.fastq.gz >> $OUTDIR/${SPE}_${IND}_R2_concat.fastq
}

### Parallelize over individuals
export INDIR OUTDIR
export -f concat
parallel --colsep '_' 'concat {1} {2}' :::: $OUTDIR/samples.txt
 