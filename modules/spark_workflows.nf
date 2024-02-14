process FASTQ_TO_SAM {
	tag "${sample}"
	label 'gatk'
	label 'mem_96GB'
	label 'core_36'
	input:
		tuple val(sample), val(is_germline), path(read_1), path(read_2)
		path(fasta)
	output:
		tuple val(sample), val(is_germline), path("${sample}_unmapped.bam"), path("${sample}_unmapped.bam.bai"), emit: ubam
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
		path(interval_list)
        path(fasta)
	output:
		tuple val(sample), val(is_germline), path("${sample}_pre_sorted.bam"), path("${sample}_pre_sorted.bam.bai"), path("${sample}_pre_sorted.bam.sbi"), emit: bam
	script:
		"""

		gatk BwaSpark \
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
		-o ${sample}_pre_sorted \
		${sample}_temp.bam

		gatk ValidateSamFile \
			-I ${sample}_pre_sorted.bam \
			-MODE SUMMARY

		gatk CreateHadoopBamSplittingIndex \
			--tmp-dir ${params.spark_tmp} \
			--create-bai \
			-I ${sample}_pre_sorted.bam \
			-O ${sample}_pre_sorted.bam.sbi 

		"""
}

process BQSR_SPARK  {
	tag "${sample}"
	label 'gatk'
	label 'mem_256GB'
	label 'core_64'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(interval_list)
		path(snv_resource)
		path(cnv_resource)
		path(sv_resource)
	output:
		tuple val(sample), path("${sample}_recal.bam"), path("${sample}_recal.bam.bai")
	script:
		"""

		ln -sf \$( echo "\$( realpath ${snv_resource} ).tbi" ) .
		ln -sf \$( echo "\$( realpath ${cnv_resource} ).tbi" ) .
		ln -sf \$( echo "\$( realpath ${sv_resource} ).tbi" ) .

		gatk CreateHadoopBamSplittingIndex \
			--create-bai \
			-I ${bam} \
			-O ${bam}.sbi 

		gatk BQSRPipelineSpark \
		--spark-master local[${task.cpus}] \
		--tmp-dir ${params.spark_tmp} \
		-L ${interval_list} \
		-I ${bam} \
		-R ${fasta}/${fasta}.fa \
		--known-sites ${snv_resource} \
		--known-sites ${cnv_resource} \
		--known-sites ${sv_resource} \
		-O ${sample}_recal.bam

		gatk ValidateSamFile \
			-I ${sample}_recal.bam \
			-MODE SUMMARY

		if \$( samtools flagstat ${sample}_recal.bam | grep -q "^0 + 0 in total (QC-passed reads + QC-failed reads)\$" ); then
        	echo "${sample}_recal.bam is empty"
			exit 1
		fi
		
		samtools index -@ ${task.cpus} ${sample}_recal.bam

		"""
}

process MARK_DUPLICATES_SPARK {
	publishDir "${params.outfolder}/${params.runID}/BAMQC", pattern: "${sample}.dup_metrics.*", mode: 'copy', overwrite: true
	publishDir "${params.outfolder}/${params.runID}/BAM", pattern: "${sample}_recal_dupfiltered.bam*", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_256GB'
	label 'core_64'
	input:
		tuple val(sample), path(bam), path(bai)
		path(interval_list)
	output:
		tuple val(sample),  path("${sample}_recal_dupfiltered.bam"), path("${sample}_recal_dupfiltered.bam.bai"), emit: ch_bam
		path("${sample}.txt"), emit: sample_checkpoint
	script:
		"""
       gatk MarkDuplicatesSpark \
			--tmp-dir ${params.spark_tmp} \
            -I ${bam} \
            -O ${sample}_recal_dupfiltered.bam \
            -M ${sample}.dup_metrics.txt \
            --conf 'spark.executor.cores=${task.cpus}'

		gatk ValidateSamFile \
			-I ${sample}_recal_dupfiltered.bam \
			-MODE SUMMARY

		if \$( samtools flagstat ${sample}_recal_dupfiltered.bam | grep -q "^0 + 0 in total (QC-passed reads + QC-failed reads)\$" ); then
        	echo "${sample}_recal_dupfiltered.bam is empty"
			exit 1
		fi

	   gatk CreateHadoopBamSplittingIndex \
			-I  ${sample}_recal_dupfiltered.bam \
			-O ${sample}_recal_dupfiltered.bam.sbi \
			--create-bai

		echo -e "${sample}\\t${params.outfolder}/${params.runID}/BAM/${sample}_recal_dupfiltered.bam\\t${params.outfolder}/${params.runID}/BAM/${sample}_recal_dupfiltered.bam.bai" > ${sample}.txt
		"""
}