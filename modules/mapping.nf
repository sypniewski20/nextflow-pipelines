process BWA_MAP_READS {
	publishDir "${params.outfolder}/${params.runID}/BAM", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_96GB'
	label 'core_24'
	input:
		tuple val(sample), path(read_1), path(read_2)
		path(fasta)
    	path(fasta_fai)
    	path(fasta_sa)
    	path(fasta_bwt)
    	path(fasta_ann)
    	path(fasta_amb)
    	path(fasta_pac)
	output:
		tuple val(sample), path("${sample}_markdups.bam"), path("${sample}_markdups.bam.bai")
	script:
		"""

		READ_GROUP_NAME=\$( zgrep @ ${read_1} | \
					head -n 1 | \
					cut -d :  -f1 | \
					sed 's/@//g' )

		tags=\$(echo "@RG\\tID:\${READ_GROUP_NAME}\\tSM:${sample}\\tPL:ILLUMINA\\tLB:lb1")

		bwa mem -t ${task.cpus} \
			-R \${tags} \
			${fasta} \
			${read_1} \
			${read_2} | \
		samblaster | \
		samtools view --reference ${fasta} \
					  --threads ${task.cpus} \
					  -b | \
		samtools sort -@ ${task.cpus} \
					  -O bam \
					  --reference ${fasta} >  ${sample}_markdups.bam

		sambamba index --nthreads ${task.cpus} ${sample}_markdups.bam

		gatk ValidateSamFile -I  ${sample}_markdups.bam \
							 -MODE	SUMMARY

		"""
}

process BWAMEM2_MAP_READS {
	publishDir "${params.outfolder}/${params.runID}/BAM", mode: 'copy', overwrite: false
	tag "${sample}"
	label 'gatk'
	label 'mem_96GB'
	label 'core_24'
	input:
		tuple val(sample), path(read_1), path(read_2)
		path(fasta)
    	path(fasta_fai)
    	path(fasta_sa)
    	path(fasta_bwt)
    	path(fasta_ann)
    	path(fasta_amb)
    	path(fasta_pac)
		path(fasta_0123)
		path(fasta_2bit)
	output:
		tuple val(sample), path("${sample}_markdups.bam"), path("${sample}_markdups.bam.bai")
	script:
		"""

		READ_GROUP_NAME=\$( zgrep @ ${read_1} | \
					head -n 1 | \
					cut -d :  -f1 | \
					sed 's/@//g' )

		tags=\$(echo "@RG\\tID:\${READ_GROUP_NAME}\\tSM:${sample}\\tPL:ILLUMINA\\tLB:lb1")

		bwa-mem2 mem -t ${task.cpus} \
			-R \${tags} \
			${fasta} \
			${read_1} \
			${read_2} | \
		samblaster | \
		samtools view --reference ${fasta} \
					  --threads ${task.cpus} \
					  -b | \
		samtools sort -@ ${task.cpus} \
					  -O bam \
					  --reference ${fasta} >  ${sample}_markdups.bam

		sambamba index --nthreads ${task.cpus} ${sample}_markdups.bam


		gatk ValidateSamFile -I  ${sample}_markdups.bam \
							 -MODE	SUMMARY

		"""
}
