process DELAMP {
    tag "${sample}"
    label 'gatk'
	label 'mem_32GB'
	label 'core_16'
	input:
		tuple val(sample), path(manta), path(delly)
	output:
		tuple val(sample), path("${sample}_cnv_paired_end_delamped.vcf.gz"), path("${sample}_cnv_paired_end_delamped.vcf.gz.tbi")
	script:
		"""

        bedtools intersect --sorted -r 0.75 -a ${manta}_sorted.vcf.gz -b ${delly}_sorted.vcf.gz | \
        bgzip -@ ${task.cpus} -c > ${sample}_cnv_paired_end_delamped.vcf.gz
        tabix -p vcf ${sample}_cnv_paired_end_delamped.vcf.gz

		"""

}

process INTERSECT_CALLS {
    tag "${sample}"
    label 'gatk'
	label 'mem_32GB'
	label 'core_16'
	input:
		tuple val(sample), path(vcf1), path(vcf2)
	output:
		tuple val(sample), path("${sample}_cnv_intersect.vcf.gz"), path(" ${sample}_cnv_intersect.vcf.gz.tbi")
	script:
		"""

        bedtools intersect --sorted -r 0.75 -a ${vcf1} -b ${vcf2} | \
        bgzip -@ ${task.cpus} -c > ${sample}_cnv_intersect.vcf.gz
        tabix -p vcf ${sample}_cnv_intersect.vcf.gz

		"""

}