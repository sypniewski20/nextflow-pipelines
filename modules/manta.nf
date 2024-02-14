process MANTA_CNV_CALL {
	tag "${sample}"
	label 'manta'
	label 'mem_4GB'
	label 'core_18'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
	output:
		tuple val(sample), path("manta/results/variants/diploidSV.vcf.gz")
	script:
		"""
            
        /manta/bin/configManta.py \
		--bam ${bam} \
		--referenceFasta ${fasta}/${fasta}.fa \
		--runDir manta

		manta/runWorkflow.py -j ${task.cpus}

		"""
}

process MANTA_FILTER_VCF {
	publishDir "${params.outfolder}/${params.runID}/CNV", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_8GB'
	label 'core_4'
	input:
		tuple val(sample), path(manta)
	output:
		tuple val(sample), path("${sample}_manta_cnv_sorted.vcf.gz"), path("${sample}_manta_cnv_sorted.vcf.gz.tbi")
	script:
		"""
            
		bcftools sort ${manta} -Ov | \
		bcftools view -f PASS --threads ${task.cpus} -Oz -o ${sample}_manta_cnv_sorted.vcf.gz

		tabix -p vcf ${sample}_manta_cnv_sorted.vcf.gz
		"""
}
