process WRITE_BAM_CHECKPOINT {
	publishDir "${params.outfolder}/${params.runID}/BAM", pattern: "bam_checkpoint.tsv", mode: 'copy', overwrite: true
	label 'wgs_tools'
	label 'mem_1GB'
	label 'core_1'
	input:
		path(samples)
	output:
		file("bam_checkpoint.tsv")
	script:
		"""

		echo -e "sampleID\\tbam\\tbai" > bam_checkpoint.tsv
		cat ${samples} >> bam_checkpoint.tsv

		"""
}