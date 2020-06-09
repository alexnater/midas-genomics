#!/bin/bash
#$ -N angsd_sfs
#$ -t 1-24:1
#$ -l h_rt=72:00:00
#$ -pe smp 4
#$ -l h_vmem=8G

function join { local IFS="$1"; shift; echo "$*"; }
export PATH=$HOME/bin:$HOME/libs/bin:$PATH

RUN_ID="${JOB_ID}.${SGE_TASK_ID}"
echo "Starting job ${JOB_ID}, task ${SGE_TASK_ID}, run ${RUN_ID}."


BIN="/data/scc3/shared/software/angsd-0.929"
WORKDIR="/data/scc3/anater/nobackup/output/angsd"
SCRATCH="/data/scc3/scratch/$RUN_ID"
REFSEQ="/data/scc3/projects/Midas_WGS/reference_genomes/Amphilophus_citrinellus.PacBio.scrubbed.marvel.BioNano.juicer.arrow.freebayes.renamed.masked.fasta"
ANCREF="/data/scc3/projects/Midas_WGS/reference_genomes/ancestral/Amphilophus_citrinellus.ancestral.mincov5_mrca.masked.fasta"
REGIONFILE="/data/scc3/anater/output/angsd/valid_regions.1kb.bed"
POPLIST="/data/scc3/anater/poplists/poplist.Midas.kautt.lips.coverage.txt"

POPS=("Nic_cit" "Nic_lab" "Apo_fla" "Apo_glo")
NIND=10
MININD_PROP=0.8
MINCOVERAGE=5
CHROMOSOME=`awk -v line=$SGE_TASK_ID 'NR==line{print $1; exit}' ${REFSEQ}.fai`
REGION="-r ${CHROMOSOME}:"
NBOOTSTRAP=100

POPSTRING=$( join _ ${POPS[@]} )
OUTFOLDER="$WORKDIR/${POPSTRING}_${NIND}ind_sfs_noexons"

BAMFOLDER="/data/scc3/projects/Midas_WGS/bam"
BAMSUFFIX=".markadap.mapped.markdup.bam"

mkdir -p $OUTFOLDER
mkdir -p $SCRATCH

cd $SCRATCH

# produce file with regions to work on:
awk -v chrom=$CHROMOSOME '$1==chrom{print $1":"$2+1"-"$3}' $REGIONFILE > regions.txt


for POP in ${POPS[@]}; do

  # create bamlist:
  awk -v pop=$POP -v path=$BAMFOLDER -v suffix=$BAMSUFFIX '$2==pop{print path "/" $1 suffix "\t" $3}' $POPLIST | sort -k2,2rn | head -n $NIND | cut -f1 > bamlist.${POP}.txt

  # figure out minimum number of individuals:
  MININD=`wc -l bamlist.${POP}.txt | awk -v minprop=$MININD_PROP '{print int($1 * minprop)}'`
  FILTERS="-minMapQ 30 -minQ 20 -minInd $MININD -minIndDepth $MINCOVERAGE"
  
  # generate the site stats:
  DOS="-doMajorMinor 1 -doMaf 1 -doSnpStat 1 -doGeno 3 -doPost 2 -doHWE 1"
  $BIN/angsd -b bamlist.${POP}.txt -gl 2 -rf regions.txt $FILTERS -skipTriallelic 1 -sb_pval 0.01 -hwe_pval 0.01 $DOS -P 4 -out ${POP}.${CHROMOSOME}.sfsSites
  
  # prepare list of filtered sites:
  zcat ${POP}.${CHROMOSOME}.sfsSites.snpStat.gz | tail -n+2 | cut -f 1,2  > ${POP}.${CHROMOSOME}.filtered_sites.txt
  sleep 60	# I have issues with the time stamp otherwise.
  $BIN/angsd sites index ${POP}.${CHROMOSOME}.filtered_sites.txt
  
  # generate the saf file, reads-in all entire chromosome, but restricts analysis to -sites:
  $BIN/angsd -b bamlist.${POP}.txt -sites ${POP}.${CHROMOSOME}.filtered_sites.txt -gl 2 -anc $ANCREF $REGION $FILTERS -doSaf 1 -P 4 -out ${POP}.${CHROMOSOME}
  
done


# estimate the 4d sfs:
$BIN/misc/realSFS ${POPS[0]}.${CHROMOSOME}.saf.idx ${POPS[1]}.${CHROMOSOME}.saf.idx ${POPS[2]}.${CHROMOSOME}.saf.idx ${POPS[3]}.${CHROMOSOME}.saf.idx -P 4 > $OUTFOLDER/${POPS[0]}_${POPS[1]}_${POPS[2]}_${POPS[3]}.${CHROMOSOME}.sfs

rm -rf $SCRATCH

