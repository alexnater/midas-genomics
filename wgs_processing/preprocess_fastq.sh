#!/bin/bash
#$ -N preprocess
#$ -t 1-482:1
#$ -tc 100
#$ -l h_rt=72:00:00
#$ -l h_vmem=16G


export TMPDIR=/data/scc3/scratch/tmp
JAVA=/usr/lib64/jvm/jre-1.8.0/bin/java

RUN_ID="${JOB_ID}.${SGE_TASK_ID}"
echo "Starting job ${JOB_ID}, task ${SGE_TASK_ID}, run ${RUN_ID}."

WORKDIR="/data/scc3/anater/output/bwa"
SCRATCH="/data/scc3/scratch/$RUN_ID"

PICARD_JAR="/data/scc3/shared/software/bin/picard.jar"

SAMPLELIST="$WORKDIR/samplelists/samplelist.all.txt"
FASTQFOLDER="/data/scc0/rawdata/Midas/WGS/fastq"
OUTFOLDER="/data/scc0/rawdata/Midas/WGS/ubam"

if [ ! -d $OUTFOLDER ]; then mkdir -p $OUTFOLDER; fi
if [ ! -d $SCRATCH ]; then mkdir -p $SCRATCH; fi

SAMPLEID=`awk -v line=$SGE_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READINFO=(`zcat $FASTQFOLDER/${SAMPLEID}_1.fq.gz | awk '/^@/{split($1,a,":"); split($2,b,":"); print a[3],a[4],b[4]; exit}'`)
ID="${READINFO[0]}.${READINFO[1]}.${READINFO[2]}"
PU="${READINFO[0]}.${READINFO[1]}.${READINFO[2]}"

$JAVA -Djava.io.tmpdir=$SCRATCH -XX:ParallelGCThreads=4 -Xmx7G -jar $PICARD_JAR FastqToSam \
	FASTQ=$FASTQFOLDER/${SAMPLEID}_1.fq.gz \
	FASTQ2=$FASTQFOLDER/${SAMPLEID}_2.fq.gz \
	OUTPUT=/dev/stdout \
	READ_GROUP_NAME=$ID \
	SAMPLE_NAME=$SAMPLEID \
	LIBRARY_NAME=$SAMPLEID \
	PLATFORM_UNIT=$PU \
	PLATFORM=illumina \
	TMP_DIR=$SCRATCH | \
	$JAVA -Djava.io.tmpdir=$SCRATCH -XX:ParallelGCThreads=4 -Xmx7G -jar $PICARD_JAR MarkIlluminaAdapters \
		INPUT=/dev/stdin \
		OUTPUT=$OUTFOLDER/${SAMPLEID}.unmapped.markadap.bam \
		METRICS=$OUTFOLDER/${SAMPLEID}.unmapped.markadap.metric \
		TMP_DIR=$SCRATCH
	

rm -rf $SCRATCH

