#!/bin/bash
#$ -N diffstats
#$ -t 1-4:1
#$ -l h_rt=24:00:00
#$ -l h_vmem=4G

function join { local IFS="$1"; shift; echo "$*"; }
export PATH=$HOME/bin:$PATH

COUNTER=$(( $SGE_TASK_ID - 1 ))
RUN_ID="${JOB_ID}.${SGE_TASK_ID}"
echo "Starting job ${JOB_ID}, task ${SGE_TASK_ID}, run ${RUN_ID}."

BIN="/data/scc3/anater/output/diffstats"
WORKDIR="/data/scc3/anater/nobackup/output/diffstats"

POPS=("Mas_cit_nlp" "Mas_cit_lip")
POPSTRING=$( join _ ${POPS[@]} )
NINDIVIDUALS=20
MIN_INDCOVERAGE=5
MAX_MISSINGNESS=0.20
UNMASKED_CUTOFF=0.2
MAX_NMISSING=2
PHASED=1

NHAPS=$((NINDIVIDUALS * 2))
WINDOWSIZE=1050000
NSUBWINDOWS=21
RHOTOTHETA=3

OUTFOLDER="$WORKDIR/data/${POPSTRING}_${NINDIVIDUALS}ind/mincov${MIN_INDCOVERAGE}_miss${MAX_MISSINGNESS}_masked"
SIMFOLDER="$OUTFOLDER/stats_msms_win${WINDOWSIZE}_nsub${NSUBWINDOWS}_rhototheta${RHOTOTHETA}_nmiss${MAX_NMISSING}"

MASKFILE="$OUTFOLDER/masks/mask.combined.fasta"
BEDFILE="none"
TRAINDIR="$WORKDIR/data/${POPSTRING}_${NINDIVIDUALS}ind/sims_msms_win${WINDOWSIZE}_rhototheta${RHOTOTHETA}/train"
TESTDIR="$WORKDIR/data/${POPSTRING}_${NINDIVIDUALS}ind/sims_msms_win${WINDOWSIZE}_rhototheta${RHOTOTHETA}/test"

mkdir -p $SIMFOLDER/train
mkdir -p $SIMFOLDER/test


cd $TRAINDIR

SIMFILES=( *.msOut )
SIMFILE=${SIMFILES[$COUNTER]}

SEED=`od -An -N2 -i /dev/random`

echo "$BIN/diffstats_masking $MASKFILE $BEDFILE $WINDOWSIZE $NSUBWINDOWS 0 $NHAPS $NHAPS $UNMASKED_CUTOFF $MAX_NMISSING $PHASED $SEED < $TRAINDIR/$SIMFILE > $SIMFOLDER/train/${SIMFILE%.msOut}.combined.txt"
$BIN/diffstats_masking $MASKFILE $BEDFILE $WINDOWSIZE $NSUBWINDOWS 0 $NHAPS $NHAPS $UNMASKED_CUTOFF $MAX_NMISSING $PHASED $SEED < $TRAINDIR/$SIMFILE > $SIMFOLDER/train/${SIMFILE%.msOut}.combined.txt

SEED=`od -An -N2 -i /dev/random`

echo "$BIN/diffstats_masking $MASKFILE $BEDFILE $WINDOWSIZE $NSUBWINDOWS 0 $NHAPS $NHAPS $UNMASKED_CUTOFF $MAX_NMISSING $PHASED $SEED < $TESTDIR/$SIMFILE > $SIMFOLDER/test/${SIMFILE%.msOut}.combined.txt"
$BIN/diffstats_masking $MASKFILE $BEDFILE $WINDOWSIZE $NSUBWINDOWS 0 $NHAPS $NHAPS $UNMASKED_CUTOFF $MAX_NMISSING $PHASED $SEED < $TESTDIR/$SIMFILE > $SIMFOLDER/test/${SIMFILE%.msOut}.combined.txt


