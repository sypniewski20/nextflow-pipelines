process GET_MT {
    tag "${sample}"
    label 'gatk'
	label 'mem_16GB'
	label 'core_16'
	input:
		tuple val(sample), path(bam), path(bai)
        path(fasta)
	output:
		tuple val(sample), path("${sample}_mt.bam"), path("${sample}_mt.bam.bai")
	script:
		"""

        samtools view ${bam} chrMT -bo ${sample}_mt.temp.bam
        samtools view -H ${sample}_mt.temp.bam | sed 's/chrMT/chrM/g' > new_header.sam
        samtools reheader new_header.sam ${sample}_mt.temp.bam > ${sample}_mt.bam
		sambamba index --nthreads ${task.cpus} ${sample}_mt.bam

		"""
}

process MT_CALL {
    publishDir "${params.outfolder}/${params.runID}/SNV/mity", mode: 'copy', overwrite: true
    label 'mity'
	label 'mem_16GB'
	label 'core_1'
	input:
		path(bam)
		path(bai)
	output:
		path("mity.vcf.gz")
	script:
		"""

        call \
        --prefix mity \
        --output-dir . \
        --normalise \
		--reference hg38 \
        ${bam}

		"""

}

process MT_REPORT {
    publishDir "${params.outfolder}/${params.runID}/SNV/mity", mode: 'copy', overwrite: true
    label 'mity'
	label 'mem_16GB'
	label 'core_1'
	input:
		path(vcf)
	output:
		tuple path("mity.annotated_variants.csv"), path("mity.annotated_variants.xlsx")
	script:
		"""

        report \
        --prefix mity \
        --min_vaf 0.01 \
        --output-dir . \
        ${vcf}

		"""

}

process MT_MERGE {
    publishDir "${params.outfolder}/${params.runID}/SNV/", mode: 'copy', overwrite: true
    label 'mity'
	label 'mem_16GB'
	label 'core_1'
	input:
		path(mity_vcf)
		path(nuclear_vcf)
	output:
		tuple path("mity.annotated_variants.csv"), path("mity.annotated_variants.xlsx")
	script:
		"""

		merge \
		--prefix multisample \
		--mity_vcf ${mity_vcf} \
		--nuclear_vcf ${nuclear_vcf}

		"""

}