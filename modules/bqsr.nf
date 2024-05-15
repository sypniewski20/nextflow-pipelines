process BASE_RECALIBRATOR {
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
    	path(fasta_fai)
		path(fasta_dict)
		path(interval_list)
		path(snv_resource)
		path(snv_tbi) 
	output:
		tuple val(sample), path("${bam}"), path("${sample}_recal_data.table")
	script:
		"""
						
		gatk BaseRecalibrator \
		-I ${bam} \
		-R ${fasta} \
		-L ${interval_list} \
		--known-sites ${snv_resource} \
		-O ${sample}_recal_data.table

		"""
}

process APPLY_BQSR {
	publishDir "${params.outfolder}/${params.runID}/BAM", mode: 'copy', overwrite: false
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_8'
	input:
		tuple val(sample), path(bam), path(recal_table)
		path(interval_list)
		path(fasta)
    	path(fasta_fai)
		path(fasta_dict)
	output:
		tuple val(sample), path("${sample}_recal.bam"), path("${sample}_recal.bam.bai")
	script:
		"""

		gatk ApplyBQSR \
		-R ${fasta} \
		-L ${interval_list} \
		-I ${bam} \
		--bqsr-recal-file ${recal_table} \
		-O ${sample}_recal.bam

		samtools index -@ ${task.cpus} ${sample}_recal.bam

		gatk ValidateSamFile -I ${sample}_recal.bam \
							 -M SUMMARY

		"""
}

process BQSR_SPARK  {
	publishDir "${params.outfolder}/${params.runID}/BAM", mode: 'copy', overwrite: false
	tag "${sample}"
	label 'gatk'
	label 'mem_64GB'
	label 'core_36'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_img)
		path(interval_list)
		path(snv_resource)
		path(snv_tbi)
	output:
		tuple val(sample), path("${sample}_sorted_markdup_recal.bam"), path("${sample}_sorted_markdup_recal.bam.bai")
	script:
		"""

		gatk CreateHadoopBamSplittingIndex \
			--create-bai \
			-I ${bam} \
			-O ${bam}.sbi 

		gatk BQSRPipelineSpark \
		--spark-master local[${task.cpus}] \
		--tmp-dir ${params.spark_tmp} \
		-L ${interval_list} \
		-I ${bam} \
		-R ${fasta} \
		--known-sites ${snv_resource} \
		-O ${sample}_sorted_markdup_recal.bam

		gatk ValidateSamFile \
			-I ${sample}_sorted_markdup_recal.bam \
			-MODE SUMMARY

		if \$( samtools flagstat ${sample}_sorted_markdup_recal.bam | grep -q "^0 + 0 in total (QC-passed reads + QC-failed reads)\$" ); then
        	echo "${sample}_recal.bam is empty"
			exit 1
		fi
		
		samtools index -@ ${task.cpus} ${sample}_recal.bam

		"""
}