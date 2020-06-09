#!/bin/bash
#$ -N run_fsc
#$ -q scc
#$ -t 1-100:1
#$ -l h_rt=48:00:00
#$ -pe smp 2
#$ -l h_vmem=4G


function join { local IFS="$1"; shift; echo "$*"; }
export PATH=$HOME/bin:$HOME/libs/bin:$PATH
export TMPDIR=/data/scc3/scratch/tmp

RUN_ID="${JOB_ID}.${SGE_TASK_ID}"
echo "Starting job ${JOB_ID}, task ${SGE_TASK_ID}, run ${RUN_ID}."

BIN="$HOME/bin/scripts"
WORKDIR="/data/scc3/anater/nobackup/output/fastsimcoal"
SCRATCH="/data/scc3/scratch/$RUN_ID"
MODEL="Man_Nic_Mas_3pops_admix"
MODELFOLDER="/data/scc3/anater/output/fastsimcoal/models/$MODEL"
TPLFILE="${MODEL}.tpl"
ESTFILE="${MODEL}.est"
NSIMS=200000
NECM=100

OUTFOLDER="${MODEL}_nofixed_noexons_intense"

mkdir -p $WORKDIR/$OUTFOLDER
mkdir -p $SCRATCH


# Attention: Does not work with absolute paths for input files!
rsync -avhSP $MODELFOLDER/ $SCRATCH

cd $SCRATCH

echo "fsc26 -t $TPLFILE -n $NSIMS -e $ESTFILE -d -M -L $NECM --multiSFS --inf -c 2 -B 2 -q"
fsc26 -t $TPLFILE -n $NSIMS -e $ESTFILE -d -M -L $NECM --multiSFS --inf -c 2 -B 2 -q

rsync -avhSP $SCRATCH/$MODEL/ $WORKDIR/$OUTFOLDER/run${SGE_TASK_ID}

rm -rf $SCRATCH

