process DELLY_CNV_CALL {
	tag "${sample}"
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
	output:
		tuple val(sample), path("${sample}_delly.vcf")
	script:
		"""

		export OMP_NUM_THREADS=${task.cpus}

        delly call -g ${fasta}/${fasta}.fa ${bam} > ${sample}_delly.vcf

		"""
}

process DELLY_FILTER_VCF {
	publishDir "${params.outfolder}/${params.runID}/CNV", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_8GB'
	label 'core_4'
	input:
		tuple val(sample), path(delly)
	output:
		tuple val(sample), path("${sample}_delly_cnv_sorted.vcf.gz"), path("${sample}_delly_cnv_sorted.vcf.gz.tbi")
	script:
		"""
            
		bcftools sort ${delly} -Ov | \
		bcftools view -f PASS --threads ${task.cpus} -Oz -o ${sample}_delly_cnv_sorted.vcf.gz

		tabix -p vcf ${sample}_delly_cnv_sorted.vcf.gz
		"""
}
