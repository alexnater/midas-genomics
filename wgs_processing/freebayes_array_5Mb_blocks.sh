#!/bin/bash
#$ -N freebayes
#$ -t 1-482:1
#$ -tc 100
#$ -l h_rt=144:00:00
#$ -l h_vmem=16G


export PATH=/data/scc3/shared/software/bin:$PATH
export TMPDIR=/data/scc3/scratch/tmp

module load gnu

RUN_ID="${JOB_ID}.${SGE_TASK_ID}"
echo "Starting job ${JOB_ID}, task ${SGE_TASK_ID}, run ${RUN_ID}."

WORKDIR="/data/scc3/anater/output/freebayes"

POP="all_Midas_withCenSiqTFA"
REGIONLIST="$WORKDIR/windows_5Mb.bed"
BAMLIST="$WORKDIR/samplelists/bamlist.all.txt"
POPLIST="$WORKDIR/poplists/poplist.all.txt"
REFSEQ="/data/scc3/projects/Midas_WGS/reference_genomes/Amphilophus_citrinellus.PacBio.scrubbed.marvel.BioNano.juicer.arrow.freebayes.renamed.fasta"
OUTFOLDER="/data/scc3/nobackup/anater/output/freebayes"

mkdir -p $OUTFOLDER

REGION=(`awk -v line=$SGE_TASK_ID 'NR==line{print $1,$2,$3; exit}' $REGIONLIST`)


freebayes -f $REFSEQ -L $BAMLIST --populations $POPLIST -r ${REGION[0]}:${REGION[1]}-${REGION[2]} -v $OUTFOLDER/${POP}.raw.snps.indels.${REGION[0]}_${REGION[1]}_${REGION[2]}.vcf --standard-filters 

