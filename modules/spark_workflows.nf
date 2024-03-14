process FASTQ_TO_SAM {
	tag "${sample}"
	label 'gatk'
	label 'mem_96GB'
	label 'core_36'
	input:
		tuple val(sample), path(read_1), path(read_2)
		path(fasta)
	output:
		tuple val(sample), path("${sample}_unmapped.bam"), path("${sample}_unmapped.bam.bai")
	script:
		"""

		READ_GROUP_NAME=\$( zgrep @ ${read_1} | \
                           head -n 1 | \
                           cut -d :  -f1 | \
                           sed 's/@//g' )

		gatk FastqToSam \
			--FASTQ ${read_1} \
			--FASTQ2 ${read_2} \
			--R ${fasta}/${fasta}.fa \
			--OUTPUT ${sample}_temp_unmapped.bam \
			--READ_GROUP_NAME \$( echo \$READ_GROUP_NAME ) \
			--SAMPLE_NAME ${sample} \
			--LIBRARY_NAME ${sample} \
			--PLATFORM_UNIT A \
			--PLATFORM illumina \
			--SEQUENCING_CENTER BI \
			--RUN_DATE ${params.date} \
			--MAX_RECORDS_IN_RAM 5000000

		sambamba view --nthreads ${task.cpus} -H ${sample}_temp_unmapped.bam > new_header.bam
		tail -n +2 ${fasta}/${fasta}.dict >> new_header.bam

		gatk ReplaceSamHeader \
			--R ${fasta}/${fasta}.fa \
			--INPUT ${sample}_temp_unmapped.bam \
			--HEADER new_header.bam \
			--OUTPUT ${sample}_unmapped.bam
		
		sambamba index --nthreads ${task.cpus} ${sample}_unmapped.bam

		"""
}
process BWA_SPARK_MAP_READS {
	tag "${sample}"
	label 'gatk'
	label 'mem_256GB'
	label 'core_64'
	input:
		tuple val(sample), path(bam), path(bai)
        path(fasta)
		path(interval_list)
	output:
		tuple val(sample), path("${sample}_sorted_markdup.bam"), path("${sample}_sorted_markdup.bam.bai"), emit: bam
	script:
		"""

		gatk BwaAndMarkDuplicatesPipelineSpark \
			--spark-master local[${task.cpus}] \
			--tmp-dir ${params.spark_tmp} \
			-R ${fasta}/${fasta}.fa \
		-L ${interval_list} \
			-I ${bam} \
			-O ${sample}_temp.bam
		
		if \$( samtools flagstat ${sample}_temp.bam | grep -q "^0 + 0 in total (QC-passed reads + QC-failed reads)\$" ); then
        	echo "${sample}_temp.bam is empty"
			exit 1
		fi

		sambamba sort --nthreads ${task.cpus} --memory-limit ${task.memory} \
		-o ${sample}_sorted_markdup \
		${sample}_temp.bam

		gatk ValidateSamFile \
			-I ${sample}_sorted_markdup.bam \
			-MODE SUMMARY

		"""
}