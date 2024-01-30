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
		--haploid_contigs="${params.haploid_contigs}" \
		--par_regions_bed=${params.par_regions_bed} \
        --reads=${bam} \
        --output_vcf=${sample}_deepvariant.vcf.gz \
        --num_shards=${task.cpus} 

		"""

}

process FILTER_SNVS {
    publishDir "${params.outfolder}/${params.runID}/SNV/", pattern: "${sample}_deepvariant_filtered", mode: 'copy', overwrite: true
    tag "${sample}"
    label 'gatk'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(vcf)
		path(fasta)
	output:
		tuple val(sample), path("${sample}_deepvariant_filtered.vcf.gz"), path("${sample}_deepvariant_filtered.vcf.gz.tbi")
	script:
		"""

		bcftools sort ${vcf} -Ov | \
		bcftools view -f PASS --threads ${task.cpus} -Oz -o ${sample}_deepvariant_filtered.vcf.gz

		tabix -p vcf ${sample}_deepvariant_filtered.vcf.gz

		"""

}