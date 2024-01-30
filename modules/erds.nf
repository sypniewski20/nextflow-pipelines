process ERDS_CNV_CALL {
	publishDir "${params.outfolder}/${params.runID}/CNV/", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'erds'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai), path(snv_calls), path(snv_calls_tbi)
		path(fasta)
	output:
		tuple val(sample), path("erds")
	script:
		"""

        erds_pipeline.pl \
        -b ${bam} \
        -v ${snv_calls} \
		-o erds \
        -r ${fasta}/${fasta}.fa

		"""
}

