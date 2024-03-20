process BWA_MAP_READS {
	tag "${sample}"
	label 'gatk'
	label 'mem_64GB'
	label 'core_36'
	input:
		tuple val(sample), path(read_1), path(read_2)
		path(fasta)
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
			${fasta}/${fasta}.fa \
			${read_1} \
			${read_2} | \
		samblaster | \
		samtools view --reference ${fasta}/${fasta}.fa \
					  --threads ${task.cpus} \
					  -b | \
		samtools sort -@ ${task.cpus} \
					  -O bam \
					  --reference ${fasta}/${fasta}.fa >  ${sample}_markdups.bam

		sambamba index --nthreads ${task.cpus} ${sample}_markdups.bam

		gatk ValidateSamFile -I  ${sample}_markdups.bam \
							 -MODE	SUMMARY

		"""
}

process BWAMEM2_MAP_READS {
	tag "${sample}"
	label 'gatk'
	label 'mem_64GB'
	label 'core_36'
	input:
		tuple val(sample), path(read_1), path(read_2)
		path(fasta)
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
			${fasta}/${fasta}.fa \
			${read_1} \
			${read_2} | \
		samblaster | \
		samtools view --reference ${fasta}/${fasta}.fa \
					  --threads ${task.cpus} \
					  -b | \
		samtools sort -@ ${task.cpus} \
					  -O bam \
					  --reference ${fasta}/${fasta}.fa >  ${sample}_markdups.bam

		sambamba index --nthreads ${task.cpus} ${sample}_markdups.bam


		gatk ValidateSamFile -I  ${sample}_markdups.bam \
							 -MODE	SUMMARY

		"""
}

process SAMBAMBA_MARK_DUPLICATES {
	publishDir "${params.outfolder}/${params.runID}/BAMQC", pattern: "${sample}.dup_metrics.*", mode: 'copy', overwrite: true
	publishDir "${params.outfolder}/${params.runID}/BAM", pattern: "${sample}_dupfiltered.*", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_8'
	input:
		tuple val(sample), path(bam), path(bai)
	output:
		tuple val(sample),  path("${sample}_dupfiltered.bam"), path("${sample}_dupfiltered.bam.bai")
	script:
		"""

		sambamba markdup \
			--nthreads ${task.cpus} \
			${bam} \
			${sample}_dupfiltered.bam
		
		sambamba index --nthreads ${task.cpus} ${sample}_dupfiltered.bam

		if \$( sambamba flagstat --nthreads ${task.cpus} ${sample}_dupfiltered.bam | grep -q "^0 + 0 in total (QC-passed reads + QC-failed reads)\$" ); then
        	echo "${sample}_dupfiltered.bam is empty"
			exit 1
		fi

		gatk ValidateSamFile -I ${sample}_dupfiltered.bam \
							 -M SUMMARY
		
		"""
}

