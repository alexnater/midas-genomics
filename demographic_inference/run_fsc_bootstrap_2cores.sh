#!/bin/bash
#$ -N fsc_boot
#$ -q scc
#$ -t 1-100:1
#$ -l h_rt=144:00:00
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
MODEL="Nic_Apo_glo_fla_4pops_admix"
MODELFOLDER="/data/scc3/anater/output/fastsimcoal/models/$MODEL"
PARFILE="${MODEL}_maxL.par"
PVFILE="${MODEL}.pv"
TPLFILE="${MODEL}.tpl"
ESTFILE="${MODEL}.est"
NREPS=10
NSIMS=200000
NECM=40

OUTFOLDER="${MODEL}_nofixed_noexons_intense"

mkdir -p $SCRATCH


# Attention: does not work with absolute paths for input files!
rsync -avhSP $MODELFOLDER/ $SCRATCH

cd $SCRATCH

# generate simulated data:
echo "fsc26 -i $PARFILE -n 1 -d -s 0 -k 10000000 --multiSFS --inf -x -c 2 -B 2 -q"
fsc26 -i $PARFILE -n 1 -d -s 0 -k 10000000 --multiSFS --inf -x -c 2 -B 2 -q

for ((i=1; i<=$NREPS; ++i)); do
  echo "Working on replicate $i ..."
  mkdir -p $SCRATCH/run${i}
  cd $SCRATCH/run${i}
  ln -s $SCRATCH/$TPLFILE $TPLFILE
  ln -s $SCRATCH/$ESTFILE $ESTFILE
  ln -s $SCRATCH/$PVFILE $PVFILE
  ln -s $SCRATCH/${PARFILE%.par}/${PARFILE%.par}_DSFS.obs ${PARFILE%_maxL.par}_DSFS.obs
  
  echo "fsc26 -t $TPLFILE -n $NSIMS -e $ESTFILE --initValues $PVFILE -d -M -L $NECM --multiSFS --inf -c 2 -B 2 -q"
  fsc26 -t $TPLFILE -n $NSIMS -e $ESTFILE --initValues $PVFILE -d -M -L $NECM --multiSFS --inf -c 2 -B 2 -q
done


awk 'NR==1{print "Run\t"$0; exit}' $SCRATCH/run1/$MODEL/${MODEL}.bestlhoods > $SCRATCH/${MODEL}.maxL.txt; for ((i=1; i<=$NREPS; ++i)); do awk -v run=$i 'NR==1{next} {print run"\t"$0}' $SCRATCH/run${i}/$MODEL/${MODEL}.bestlhoods; done >> $SCRATCH/${MODEL}.maxL.txt

rsync -avhSP $SCRATCH/${MODEL}.maxL.txt $WORKDIR/$OUTFOLDER/${MODEL}.bootstrap.${SGE_TASK_ID}.txt

rm -rf $SCRATCH

