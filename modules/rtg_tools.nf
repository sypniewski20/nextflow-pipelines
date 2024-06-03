process RTG_SDF {
    label 'rtg'
	label 'mem_2GB'
	label 'core_1'
	input:
		path(fasta)
	output:
		path("${fasta}.sdf")
	script:
		"""

        rtg format \
        -o ${fasta}.sdf \
        ${fasta}

		"""

}

process RTG_MENDELIAN {
    publishDir "${params.outfolder}/${params.runID}/SNV/", mode: 'copy', overwrite: true
    label 'rtg'
	label 'mem_2GB'
	label 'core_1'
	input:
		path(vcf)
	output:
		tuple path("annotated_${vcf}"), path("annotated_${vcf}.tbi")
	script:
		"""

        mendelian \
        -i ${vcf} \
        -o annotated_${vcf} \
        --pedigree=${params.ped} \
        -t /data/references/fasta/GRCh38/GRCh38.sdf

		index annotated_${vcf}

		"""

}