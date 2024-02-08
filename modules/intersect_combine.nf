process INTERSECT_PAIRED_END_CNV {
	publishDir "${params.outfolder}/${params.runID}/CNV", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'	
	label 'mem_8GB'
	label 'core_4'
	input:
		tuple val(sample), path(delly), path(delly_tbi), path(manta), path(manta_tbi)
	output:
		tuple val(sample), path("${sample}_pairend.vcf.gz"), path("${sample}_pairend.vcf.gz.tbi")
	script:
		"""

        bedtools intersect -f 0.75 -r -a ${delly} -b ${manta} > cnv.vcf
		bcftools sort cnv.vcf -Oz -o ${sample}_pairend.vcf.gz

		tabix -p vcf ${sample}_pairend.vcf.gz
		"""
}

process INTERSECT_COVERAGE_CNV {
	publishDir "${params.outfolder}/${params.runID}/CNV", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'	
	label 'mem_8GB'
	label 'core_4'
	input:
		tuple val(sample), path(cnvpytor), path(cnvpytor_tbi), path(erds), path(erds_tbi)
	output:
		tuple val(sample), path("${sample}_coverage.vcf.gz"), path("${sample}_coverage.vcf.gz.tbi")
	script:
		"""

        bedtools intersect -f 0.75 -r -a ${cnvpytor} -b ${erds} > cnv.vcf   
		bcftools sort cnv.vcf -Oz -o ${sample}_coverage.vcf.gz

		tabix -p vcf ${sample}_coverage.vcf.gz
		"""
}

process COMBINE_CNV {
	publishDir "${params.outfolder}/${params.runID}/CNV", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'	
	label 'mem_8GB'
	label 'core_4'
	input:
		file(vcf)
		file(tbi)
	output:
		tuple val(sample), path("multisample_cnv_combined.vcf.gz"), path("multisample_cnv_combined.vcf.gz.tbi")
	script:
		"""

		bcftools merge ${vcf} -Ou | \
		bcftools sort -Oz -o multisample_cnv_combined.vcf.gz

		tabix -p vcf multisample_cnv_combined.vcf.gz

		"""
}