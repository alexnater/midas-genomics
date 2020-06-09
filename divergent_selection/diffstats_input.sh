#!/bin/bash
#$ -N diffstats
#$ -q scc
#$ -l h_rt=24:00:00
#$ -l h_vmem=4G


function join { local IFS="$1"; shift; echo "$*"; }
export PATH=$HOME/bin:$PATH

echo "Starting job ${JOB_ID}."

BIN="/data/scc3/anater/output/diffstats"
WORKDIR="/data/scc3/anater/nobackup/output/diffstats"

POPS=("Xil_vir" "Xil_xil")
POPSTRING=$( join _ ${POPS[@]} )
NINDIVIDUALS=15
MIN_INDCOVERAGE=5
MAX_MISSINGNESS=0.20
UNMASKED_CUTOFF=0.2
MAX_NMISSING=2
PHASED=1

NHAPS=$((NINDIVIDUALS * 2))
WINDOWSIZE=1050000
NSUBWINDOWS=21

OUTFOLDER="$WORKDIR/data/${POPSTRING}_${NINDIVIDUALS}ind/mincov${MIN_INDCOVERAGE}_miss${MAX_MISSINGNESS}_masked"
INFOLDER="$OUTFOLDER/input_win${WINDOWSIZE}_nsub${NSUBWINDOWS}"

MASKFILE="$OUTFOLDER/masks/mask.combined.fasta"
INFILE="$INFOLDER/input_${POPSTRING}_${NINDIVIDUALS}ind.all.mod.txt"
OUTFILE="$INFOLDER/input_${POPSTRING}_${NINDIVIDUALS}ind.stats.nmiss${MAX_NMISSING}.txt"
BEDFILE="$INFOLDER/regions.all.bed"

SEED=`od -An -N2 -i /dev/random`

echo "$BIN/diffstats_masking $MASKFILE $BEDFILE $WINDOWSIZE $NSUBWINDOWS 0 $NHAPS $NHAPS $UNMASKED_CUTOFF $MAX_NMISSING $PHASED $SEED < $INFILE > $OUTFILE"
$BIN/diffstats_masking $MASKFILE $BEDFILE $WINDOWSIZE $NSUBWINDOWS 0 $NHAPS $NHAPS $UNMASKED_CUTOFF $MAX_NMISSING $PHASED $SEED < $INFILE > $OUTFILE

grep -v 'nan' $OUTFILE > ${OUTFILE%.txt}.cleaned.txt

