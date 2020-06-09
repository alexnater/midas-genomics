#!/bin/bash
#$ -N phase_PIRs
#$ -t 1-24:1
#$ -l h_rt=300:00:00
#$ -pe smp 4
#$ -l h_vmem=8G


export PATH=/data/scc3/shared/software/bin:$PATH
export TMPDIR=/data/scc3/scratch/tmp

module load vcftools

RUN_ID="${JOB_ID}.${SGE_TASK_ID}"
echo "Starting job ${JOB_ID}, task ${SGE_TASK_ID}, run ${RUN_ID}."

SPECIES="Midas"
WORKDIR="/data/scc3/anater/output/shapeit"
SCRATCH="/data/scc3/scratch/$RUN_ID"
REFSEQ="/data/scc3/projects/Midas_WGS/reference_genomes/Amphilophus_citrinellus.PacBio.scrubbed.marvel.BioNano.juicer.arrow.freebayes.renamed.fasta"
BAMLIST="bamlist_${SPECIES}.txt"
VCFFOLDER="/data/scc3/projects/Midas_WGS/vcf"
VCFFILE="all_Midas_withCenSiqTFA.sitefilt.snps.indels.norm.decomposed.vcf.gz"
PIRFILE_PREFIX="PIRlist"
OUTFOLDER="$WORKDIR/phasing_${SPECIES}"
OUTFILE_PREFIX="phasing"
SEED=`od -An -N3 -i /dev/random`


mkdir -p $OUTFOLDER
mkdir -p $SCRATCH

CHROMOSOME=`awk -v line=$SGE_TASK_ID 'NR==line {print $1; exit}' ${REFSEQ}.fai`
PIRFILE="${PIRFILE_PREFIX}_${SPECIES}_${CHROMOSOME}.txt"
OUTFILE="${OUTFILE_PREFIX}_${SPECIES}_${CHROMOSOME}"

if [ ! -f $SCRATCH/${BAMLIST%.txt}_${CHROMOSOME}.txt ]; then
  echo "sed -e s/"{chromosome}"/${CHROMOSOME}/ $WORKDIR/$BAMLIST > $SCRATCH/${BAMLIST%.txt}_${CHROMOSOME}.txt"
  sed -e s/"{chromosome}"/${CHROMOSOME}/ $WORKDIR/$BAMLIST > $SCRATCH/${BAMLIST%.txt}_${CHROMOSOME}.txt
fi

if [ ! -f $SCRATCH/samplelist_${SPECIES}.txt ]; then
  echo "cut -f1 $WORKDIR/$BAMLIST > $SCRATCH/samplelist_${SPECIES}.txt"
  cut -f1 $WORKDIR/$BAMLIST > $SCRATCH/samplelist_${SPECIES}.txt
fi


if [ ! -f $SCRATCH/${VCFFILE%.vcf.gz}_${CHROMOSOME}.recode.vcf ]; then
  echo "vcftools --gzvcf $VCFFOLDER/$VCFFILE --chr $CHROMOSOME --remove-indels --keep $SCRATCH/samplelist_${SPECIES}.txt --max-alleles 2 --max-missing 0.00001 --recode --out $SCRATCH/${VCFFILE%.vcf.gz}_${CHROMOSOME}"
  vcftools --gzvcf $VCFFOLDER/$VCFFILE --chr $CHROMOSOME --remove-indels --keep $SCRATCH/samplelist_${SPECIES}.txt --max-alleles 2 --max-missing 0.00001 --recode --out $SCRATCH/${VCFFILE%.vcf.gz}_${CHROMOSOME}
fi


cd $SCRATCH

if [ ! -f $SCRATCH/$PIRFILE ]; then
  echo "extractPIRs --base-quality 13 --read-quality 10 --bam ${BAMLIST%.txt}_${CHROMOSOME}.txt --vcf ${VCFFILE%.vcf.gz}_${CHROMOSOME}.recode.vcf --out $PIRFILE"
  extractPIRs --base-quality 13 --read-quality 10 --bam ${BAMLIST%.txt}_${CHROMOSOME}.txt --vcf ${VCFFILE%.vcf.gz}_${CHROMOSOME}.recode.vcf --out $PIRFILE
fi

if [ ! -f $OUTFOLDER/$OUTFILE.vcf.gz ]; then
  echo "shapeit -assemble --thread 4 --seed $SEED --states 200 --window 0.5 --burn 10 --prune 10 --main 50 --force --input-vcf ${VCFFILE%.vcf.gz}_${CHROMOSOME}.recode.vcf --input-pir $PIRFILE -O $OUTFILE"
  shapeit -assemble --thread 4 --seed $SEED --states 200 --window 0.5 --burn 10 --prune 10 --main 50 --force --input-vcf ${VCFFILE%.vcf.gz}_${CHROMOSOME}.recode.vcf --input-pir $PIRFILE -O $OUTFILE

  echo "shapeit -convert --input-haps $OUTFILE --output-vcf ${OUTFILE}.vcf"
  shapeit -convert --input-haps $OUTFILE --output-vcf ${OUTFILE}.vcf

  echo "bgzip ${OUTFILE}.vcf && tabix -p vcf ${OUTFILE}.vcf.gz"
  bgzip ${OUTFILE}.vcf && tabix -p vcf ${OUTFILE}.vcf.gz

  mv $SCRATCH/${OUTFILE}.vcf* $OUTFOLDER
fi

# to safeguard against deleting the entire scratch folder if RUN_ID is not set:
if [ ! -z $RUN_ID ]; then rm -rf $SCRATCH; fi

