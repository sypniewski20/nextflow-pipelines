process SV2_REGENOTYPE {
	publishDir "${params.outfolder}/${params.runID}/CNV/sv2", mode: 'copy', overwrite: true
	tag "${sampleID}"
	label 'sv2'	
	label 'mem_8GB'
	label 'core_4'
	input:
		tuple val(sampleID), path(bam), path(bai)
		tuple val(sampleID), path(cnv)
		tuple val(sampleID), path(snv)
		path(ped_file)
		path(fasta)
	output:
        path("*")
	script:
		"""

		for vcf in ${cnv} ${snv}; do
			ln -s "\$( realpath \${vcf} ).tbi" .
		done

        sv2 -hg38 ${fasta}/${fasta}.fa
			-i ${bam} \
			-v ${cnv} \
			-snv ${snv} \
			-p ${ped_file} \
			-o ${sampleID}_sv2

		"""
}