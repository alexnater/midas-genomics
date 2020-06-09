#!/bin/bash
#$ -N doc
#$ -l h_rt=300:00:00
#$ -pe smp 8
#$ -l h_vmem=8G


# Job submission script for Grid Engine.
# Takes mapped, merged, and duplicate-marked bam files and outputs individiual seuqencing coverage at all sites in the reference genome. Output is bgzipped and tabix-indexed.
# Requires Samtools, generate_regions.py from Freebayes package, GNU Parallel, bgzip, and tabix.


export PATH=/data/scc3/shared/software/bin:$PATH
export TMPDIR=/data/scc3/scratch/tmp

WORKDIR="/data/scc3/anater/depth"
SCRATCH="/data/scc3/scratch/$JOB_ID"

SAMPLELIST="$WORKDIR/poplists/poplist.all.txt"
REFSEQ="/data/scc3/projects/Midas_WGS/reference_genomes/Amphilophus_citrinellus.PacBio.scrubbed.marvel.BioNano.juicer.arrow.freebayes.renamed.fasta"
BAMFOLDER="/data/scc3/projects/Midas_WGS/bam"
OUTFOLDER="/data/scc3/projects/Midas_WGS/depth"
OUTFILE="samtools_depth.all.txt"
MIN_BASEQUAL=20
MIN_MAPPINGQUAL=30
REGION_SIZE=1000000

mkdir -p $OUTFOLDER
mkdir -p $SCRATCH

function join_by { local IFS="$1"; shift; echo "$*"; }

SAMPLES=(); BAMIDS=()
while read ID other; do
  SAMPLES+=($ID)
  BAMIDS+=($BAMFOLDER/${ID}.markadap.mapped.markdup.bam)
done < $SAMPLELIST


join_by $'\t' "#CHROM" "POS" ${SAMPLES[@]} > $OUTFOLDER/$OUTFILE

generate_regions.py $REFSEQ $REGION_SIZE | parallel --tmpdir $SCRATCH -k -j 16 samtools depth -q $MIN_BASEQUAL -Q $MIN_MAPPINGQUAL -r {} ${BAMIDS[@]} >> $OUTFOLDER/$OUTFILE

bgzip $OUTFOLDER/$OUTFILE && tabix -s 1 -b 2 -e 2 $OUTFOLDER/${OUTFILE}.gz

rm -rf $SCRATCH

