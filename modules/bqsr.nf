process BASE_RECALIBRATOR {
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(interval_list)
		tuple path(snv_resource), path("${snv_resource}.tbi") 
	output:
		tuple val(sample), path("${bam}"), path("${sample}_recal_data.table")
	script:
		"""
			
		ln -sf \$( echo "\$( realpath ${snv_resource} ).tbi" ) .
			
		gatk BaseRecalibrator \
		-I ${bam} \
		-R ${fasta}/${fasta}.fa \
		-L ${interval_list} \
		--known-sites ${snv_resource} \
		-O ${sample}_recal_data.table
		"""
}

process APPLY_BQSR {
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_8'
	input:
		tuple val(sample), path(bam), path(recal_table)
		path(interval_list)
		path(fasta)
	output:
		tuple val(sample), path("${sample}_recal.bam"), path("${sample}_recal.bam.bai")
	script:
		"""

		gatk ApplyBQSR \
		-R ${fasta}/${fasta}.fa \
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
	tag "${sample}"
	label 'gatk'
	label 'mem_256GB'
	label 'core_64'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(interval_list)
		path(snv_resource)
	output:
		tuple val(sample), path("${sample}_sorted_markdup_recal.bam"), path("${sample}_sorted_markdup_recal.bam.bai")
	script:
		"""

		ln -sf \$( echo "\$( realpath ${snv_resource} ).tbi" ) .

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