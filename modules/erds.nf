process ERDS_CNV_CALL {
	publishDir "${params.outfolder}/${params.runID}/CNV/", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'erds'
	label 'mem_16GB'
	label 'core_16'
	input:
		tuple val(sample), path(bam), path(bai), path(snv_calls)
		path(fasta)
	output:
		tuple val(sample), path("${sample}.erds.vcf")
	script:
		"""

        perl /ERDS/erds1.1/erds_pipeline.pl \
        -b ${bam} \
        -v ${snv_calls} \
		-o . \
        -r ${fasta}/${fasta}.fa

		"""
}

process ERDS_FILTER_VCF {
	publishDir "${params.outfolder}/${params.runID}/CNV", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'	
	label 'mem_8GB'
	label 'core_4'
	input:
		tuple val(sample), path(vcf)
	output:
		tuple val(sample), path("${sample}_erds_sorted.vcf.gz"), emit: vcf
		tuple val(sample), path("${sample}_erds_sorted.vcf.gz.tbi"), emit: tbi
	script:
		"""
            
		bcftools sort ${vcf} -Ov | \
		bcftools view -f PASS --threads ${task.cpus} -Oz -o ${sample}_erds_sorted.vcf.gz

		tabix -p vcf ${sample}_erds_sorted.vcf.gz
		"""
}

