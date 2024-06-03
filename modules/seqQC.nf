process UBAM {
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_4'
	input:
		tuple val(sample), path(bam), path(bai)
	output:
		tuple val(sample), path("${sample}_ubam.bam")
	script:
		"""

		gatk RevertSam \
		I=${bam} \
		O=${sample}_ubam.bam \
		SANITIZE=true \
		ATTRIBUTE_TO_CLEAR=XT \
		ATTRIBUTE_TO_CLEAR=XN \
		ATTRIBUTE_TO_CLEAR=AS \
		ATTRIBUTE_TO_CLEAR=OC \
		ATTRIBUTE_TO_CLEAR=OP \
		SORT_ORDER=queryname \
		RESTORE_ORIGINAL_QUALITIES=true \
		REMOVE_DUPLICATE_INFORMATION=true \
		REMOVE_ALIGNMENT_INFORMATION=true \
		MAX_DISCARD_FRACTION=0.15

		"""
}

process BAM2FASTQ {
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_4'
	input:
		tuple val(sample), path(bam)
	output:
		tuple val(sample), path("${sample}_1.fq.gz"), path("${sample}_2.fq.gz")
	script:
		"""

		gatk SamToFastq \
			NON_PF=true \
			I=${bam} \
			FASTQ=${sample}_1.fq.gz \
			F2=${sample}_2.fq.gz
 

		"""
}	

process FASTQC_PROCESSING {
	publishDir "${params.outfolder}/${params.runID}/fastQC", pattern: "fastqc_${sample}_logs/*", mode: 'copy', overwrite: true
	label 'gatk'
	tag "${sample}"
	label 'mem_8GB'
	label 'core_4'
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
	publishDir "${params.outfolder}/${params.runID}/BAMQC", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_4GB'
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