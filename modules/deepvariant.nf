process DEEP_VARIANT {
    tag "${sample}"
    label 'deepvariant'
	label 'mem_96GB'
	label 'core_32'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
	output:
		tuple val(sample), path("${sample}_deepvariant.vcf.gz")
	script:
		"""

        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=${fasta}/${fasta}.fa \
        --reads=${bam} \
        --output_vcf=${sample}_deepvariant.vcf.gz \
        --num_shards=${task.cpus} 

		"""

}

process FILTER_AND_MERGE_SNVS {
    publishDir "${params.outfolder}/${params.runID}/SNV/", mode: 'copy', overwrite: true
    tag "${sample}"
    label 'gatk'
	label 'mem_16GB'
	label 'core_8'
	input:
		tuple val(sample), path(vcf)
		path(fasta)
	output:
		tuple val(sample), path("${sample}_deepvariant_filtered.vcf.gz"), path("${sample}_deepvariant_filtered.vcf.gz.tbi"), emit: vcf 
	script:
		"""

		bcftools sort ${vcf} -Ou | \
		bcftools view -f PASS -Oz -o ${sample}_deepvariant_filtered.vcf.gz

		tabix -p vcf ${sample}_deepvariant_filtered.vcf.gz

		"""

}