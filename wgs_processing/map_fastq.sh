#!/bin/bash
#$ -N mapping
#$ -t 1-482:1
#$ -tc 100
#$ -l h_rt=72:00:00
#$ -pe smp 4
#$ -l h_vmem=8G


export PATH=/data/scc3/shared/software/bin:$PATH
export TMPDIR=/data/scc3/scratch/tmp
JAVA=/usr/lib64/jvm/jre-1.8.0/bin/java

RUN_ID="${JOB_ID}.${SGE_TASK_ID}"
echo "Starting job ${JOB_ID}, task ${SGE_TASK_ID}, run ${RUN_ID}."

WORKDIR="/data/scc3/anater/output/bwa"
SCRATCH="/data/scc3/scratch/$RUN_ID"

PICARD_JAR="/data/scc3/shared/software/bin/picard.jar"

SAMPLELIST="$WORKDIR/samplelist.all.txt"
REFSEQ="/data/scc3/projects/Midas_WGS/reference_genomes/Amphilophus_citrinellus.PacBio.scrubbed.marvel.BioNano.juicer.arrow.freebayes.renamed.fasta"
UBAMFOLDER="/data/scc0/rawdata/Midas/WGS/ubam"
OUTFOLDER="/data/scc3/projects/Midas_WGS/bam"

mkdir -p $OUTFOLDER
mkdir -p $SCRATCH
mkdir -p $SCRATCH/tmp

SAMPLEID=`awk -v line=$SGE_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`

# copy ubam file to scratch:
rsync -avhSP $UBAMFOLDER/${SAMPLEID}.unmapped.markadap.bam $SCRATCH

$JAVA -Djava.io.tmpdir=$SCRATCH/tmp -XX:ParallelGCThreads=4 -Xmx8G -jar $PICARD_JAR SamToFastq \
	I=$SCRATCH/${SAMPLEID}.unmapped.markadap.bam \
	FASTQ=/dev/stdout \
	CLIPPING_ATTRIBUTE=XT \
	CLIPPING_ACTION=2 \
	INTERLEAVE=true \
	NON_PF=true \
	TMP_DIR=$SCRATCH/tmp | \
	bwa mem -K 100000000 -t 4 -p $REFSEQ \
		/dev/stdin | \
		$JAVA -Djava.io.tmpdir=$SCRATCH/tmp -XX:ParallelGCThreads=4 -Xmx16G -jar $PICARD_JAR MergeBamAlignment \
			R=$REFSEQ \
			UNMAPPED_BAM=$SCRATCH/${SAMPLEID}.unmapped.markadap.bam \
			ALIGNED_BAM=/dev/stdin \
			OUTPUT=$SCRATCH/${SAMPLEID}.markadap.mapped.bam \
			CREATE_INDEX=false \
			ADD_MATE_CIGAR=true \
			CLIP_ADAPTERS=false \
			CLIP_OVERLAPPING_READS=true \
			INCLUDE_SECONDARY_ALIGNMENTS=true \
			MAX_INSERTIONS_OR_DELETIONS=-1 \
			PRIMARY_ALIGNMENT_STRATEGY=BestMapq\
			ATTRIBUTES_TO_RETAIN=XS \
			ATTRIBUTES_TO_RETAIN=XT \
			TMP_DIR=$SCRATCH/tmp

$JAVA -Djava.io.tmpdir=$SCRATCH/tmp -XX:ParallelGCThreads=4 -Xmx24G -jar $PICARD_JAR MarkDuplicates \
	OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
	I=$SCRATCH/${SAMPLEID}.markadap.mapped.bam \
	O=$OUTFOLDER/${SAMPLEID}.markadap.mapped.markdup.bam \
	M=$OUTFOLDER/${SAMPLEID}.markadap.mapped.markdup.metrics \
	TMP_DIR=$SCRATCH/tmp

samtools index $OUTFOLDER/${SAMPLEID}.markadap.mapped.markdup.bam

rm -rf $SCRATCH

