#!/bin/bash
#$ -N train_predict
#$ -q scc
#$ -l h_rt=6:00:00
#$ -pe smp 2
#$ -l h_vmem=4G

function join { local IFS="$1"; shift; echo "$*"; }
export PATH=/home/scc/anater/bin:/data/scc3/shared/software/anaconda3/bin:$PATH

echo "Starting job ${JOB_ID}."

BINDIR="/data/scc3/anater/output/diffstats/msms_mom"
WORKDIR="/data/scc3/anater/nobackup/output/diffstats"

POPS=("Mas_cit_nlp" "Mas_cit_lip")
POPSTRING=$( join _ ${POPS[@]} )
NINDIVIDUALS=20
MIN_INDCOVERAGE=5
MAX_MISSINGNESS=0.20
MAX_NMISSING=2

WINDOWSIZE=1050000
NSUBWINDOWS=21
RHOTOTHETA=3

KEEP_TRAIN=4000
KEEP_TEST=2000
EPOCHS=100
PATIENCE=10

MODEL_PREFIX="e${EPOCHS}_p${PATIENCE}"
OUTFOLDER="$WORKDIR/data/${POPSTRING}_${NINDIVIDUALS}ind/mincov${MIN_INDCOVERAGE}_miss${MAX_MISSINGNESS}_masked"
STATSFOLDER="$OUTFOLDER/stats_msms_win${WINDOWSIZE}_nsub${NSUBWINDOWS}_rhototheta${RHOTOTHETA}_nmiss${MAX_NMISSING}"
INPUTFOLDER="$OUTFOLDER/input_win${WINDOWSIZE}_nsub21"
INPUTFILE="$INPUTFOLDER/input_${POPSTRING}_${NINDIVIDUALS}ind.stats.nmiss${MAX_NMISSING}.cleaned.txt"
OUTFILE="$INPUTFOLDER/predictions/input_${POPSTRING}_${NINDIVIDUALS}ind.4cats.preds.${MODEL_PREFIX}.txt"

mkdir -p $STATSFOLDER/models/
mkdir -p $INPUTFOLDER/predictions/

if [ -f $STATSFOLDER/train/neutral.combined.txt ]; then
  for file in $STATSFOLDER/train/*.combined.txt; do grep -v 'nan' $file | cut -f5- > ${file%.combined.txt}.fvec; done
  rm -f $STATSFOLDER/train/*.combined.txt
fi

if [ -f $STATSFOLDER/test/neutral.combined.txt ]; then  
  for file in $STATSFOLDER/test/*.combined.txt; do grep -v 'nan' $file | cut -f5- > ${file%.combined.txt}.fvec; done
  rm -f $STATSFOLDER/test/*.combined.txt
fi


echo "python $BINDIR/train_model_4cats.py --equal --insist --keep_train $KEEP_TRAIN --keep_test $KEEP_TEST --epochs $EPOCHS --patience $PATIENCE --numSubWins $NSUBWINDOWS --predictFile $STATSFOLDER/models/${MODEL_PREFIX}.valpreds.txt $STATSFOLDER/train/ $STATSFOLDER/test/ $STATSFOLDER/models/$MODEL_PREFIX"
python $BINDIR/train_model_4cats.py --equal --insist --keep_train $KEEP_TRAIN --keep_test $KEEP_TEST --epochs $EPOCHS --patience $PATIENCE --numSubWins $NSUBWINDOWS --predictFile $STATSFOLDER/models/${MODEL_PREFIX}.valpreds.txt $STATSFOLDER/train/ $STATSFOLDER/test/ $STATSFOLDER/models/$MODEL_PREFIX

echo "python $BINDIR/predict_4cats.py --numSubWins $NSUBWINDOWS $STATSFOLDER/models/${MODEL_PREFIX}.json $STATSFOLDER/models/${MODEL_PREFIX}.weights.hdf5 $INPUTFILE $OUTFILE"
python $BINDIR/predict_4cats.py --numSubWins $NSUBWINDOWS $STATSFOLDER/models/${MODEL_PREFIX}.json $STATSFOLDER/models/${MODEL_PREFIX}.weights.hdf5 $INPUTFILE $OUTFILE

