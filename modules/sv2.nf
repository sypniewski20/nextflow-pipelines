process SV2_REGENOTYPE {
	publishDir "${params.outfolder}/${params.runID}/CNV", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'	
	label 'mem_8GB'
	label 'core_4'
	input:
		tuple val(sample), path(bam), path(cnv), path(snv), path(snv_tbi)
	output:
        path("*")
	script:
		"""

        sv2 -i ${bam} -v ${cnv} -snv ${snv} -p /home/mateuszsypniewski/S8436/S8436.ped

		"""
}