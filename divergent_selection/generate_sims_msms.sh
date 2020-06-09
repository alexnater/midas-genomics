#!/bin/bash
#$ -N generate_sims
#$ -t 1-22:1
#$ -l h_rt=48:00:00
#$ -pe smp 2
#$ -l h_vmem=4G

function join { local IFS="$1"; shift; echo "$*"; }
export PATH=/home/scc/anater/bin:/data/scc3/shared/software/anaconda3/bin:$PATH

RUN_ID="${JOB_ID}.${SGE_TASK_ID}"
echo "Starting job ${JOB_ID}, task ${SGE_TASK_ID}."

BINDIR="/data/scc3/anater/output/diffstats/msms_mom"
WORKDIR="/data/scc3/nobackup/anater/output/diffstats"

POPS=("Mas_cit_nlp" "Mas_cit_lip")
PARAMFILE="/data/scc3/anater/output/diffstats/msms_mom/models/Man_Nic_Mas_4pops_admix.maxL.txt"
POPLABELS="nic,man,nlp,lip"
POPSTRING=$( join _ ${POPS[@]} )
INDIVIDUALS=20
NHAPS=$((INDIVIDUALS * 2))

WINDOWSIZE=1050000
RHOTOTHETA=3
NSIMS_TRAIN=200
NSIMS_TEST=100

OUTFOLDER="$WORKDIR/data/${POPSTRING}_${INDIVIDUALS}ind/sims_msms_win${WINDOWSIZE}_rhototheta${RHOTOTHETA}"
OUTFILE_PREFIX="sims.${SGE_TASK_ID}"
TRAINDIR="$OUTFOLDER/train"
TESTDIR="$OUTFOLDER/test"

mkdir -p $TRAINDIR
mkdir -p $TESTDIR


cd $WORKDIR

echo "python $BINDIR/generateSimLaunchScript_msms.py $TRAINDIR $TESTDIR $OUTFILE_PREFIX $NSIMS_TRAIN $NSIMS_TEST 0,0,${NHAPS},${NHAPS} $WINDOWSIZE $RHOTOTHETA $PARAMFILE $POPLABELS > $OUTFOLDER/simfile.${SGE_TASK_ID}.sh"
python $BINDIR/generateSimLaunchScript_msms.py $TRAINDIR $TESTDIR $OUTFILE_PREFIX $NSIMS_TRAIN $NSIMS_TEST 0,0,${NHAPS},${NHAPS} $WINDOWSIZE $RHOTOTHETA $PARAMFILE $POPLABELS > $OUTFOLDER/simfile.${SGE_TASK_ID}.sh

echo "bash $OUTFOLDER/simfile.${SGE_TASK_ID}.sh"
bash $OUTFOLDER/simfile.${SGE_TASK_ID}.sh


