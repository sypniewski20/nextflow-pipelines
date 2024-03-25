process CHECK_INTEGRITY {
	label 'gatk'
	tag "${sample}"
	label 'mem_16GB'
	label 'core_8'
	input:
		tuple val(sample),  file(read_1), file(read_2)
	output:
		tuple val(sample),  file("${read_1}"), file("${read_2}")
	script:
	"""

	bgzip -@ ${task.cpus} -t ${read_1}
	bgzip -@ ${task.cpus} -t ${read_2}

	"""
}

process FASTQC_PROCESSING {
	publishDir "${params.outfolder}/${params.runID}/fastQC", pattern: "fastqc_${sample}_logs/*", mode: 'copy', overwrite: true
	label 'gatk'
	tag "${sample}"
	label 'mem_16GB'
	label 'core_8'
	input:
		tuple val(sample), file(read_1), file(read_2)
	output:
		path("fastqc_${sample}_logs/*")
	script:
		"""
		mkdir fastqc_${sample}_logs
		fastqc -o fastqc_${sample}_logs -t ${task.cpus} ${read_1} ${read_2} 
		"""
}

process FASTP_PROCESSING {
	publishDir "${params.outfolder}/${params.runID}/fastp/${sample}", pattern: "fastp.*", mode: 'copy', overwrite: true
	label 'gatk'
	tag "${sample}"
	label 'mem_8GB'
	label 'core_8'
	input:
		tuple val(sample), file(read_1), file(read_2)
	output:
		tuple val(sample), path("${sample}.filtered.R1.fq.gz"), path("${sample}.filtered.R2.fq.gz"), emit: fastq_filtered
		path("fastp.*"), emit: fastp_log
	script:
		"""
		fastp -i ${read_1} \
			  -I ${read_2} \
			  -o ${sample}.filtered.R1.fq.gz \
			  -O ${sample}.filtered.R2.fq.gz \
			  -t ${task.cpus}
		"""
}

process MOSDEPTH_WGS {
	publishDir "${params.outfolder}/${params.runID}/BAMQC", pattern: "*mosdepth*", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_4'
	input:
		tuple val(sample), path(bam), path(bai)
	output:
		tuple val(sample), path("${sample}*")
	script:
		"""

		mosdepth -n --fast-mode --by 500  ${sample} ${bam} --threshold 1,10,20,30 -t ${task.cpus}

		"""
}

process MOSDEPTH_EXOME {
	publishDir "${params.outfolder}/${params.runID}/BAMQC", pattern: "*mosdepth*", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_4'
	input:
		tuple val(sample), path(bam), path(bai)
	output:
		tuple val(sample), path("${sample}*")
	script:
		"""

		mosdepth --by ${params.contigs_bed} -n --fast-mode ${sample} ${bam} --threshold 1,10,20,30 -t ${task.cpus}

		"""
}