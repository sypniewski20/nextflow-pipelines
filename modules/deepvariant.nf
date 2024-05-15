fasta = file(params.fasta)

process DEEP_VARIANT_WGS {
    tag "${sample}"
    label 'deepvariant'
	label 'mem_96GB'
	label 'core_32'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
	output:
		tuple val(sample), path("${sample}_deepvariant.vcf.gz")
	script:
		"""

        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=${fasta}/${fasta}.fa \
        --reads=${bam} \
		--regions=${params.contigs_bed} \
        --output_vcf=${sample}_deepvariant.vcf.gz \
        --num_shards=${task.cpus} 

		"""

}

process DEEP_VARIANT_WES {
    tag "${sample}"
    label 'deepvariant'
	label 'mem_96GB'
	label 'core_32'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
	output:
		tuple val(sample), path("${sample}_deepvariant.vcf.gz")
	script:
		"""

        /opt/deepvariant/bin/run_deepvariant \
        --model_type WES \
        --ref ${fasta} \
        --reads ${bam} \
        --output_vcf ${sample}_deepvariant.vcf.gz \
        --num_shards ${task.cpus} \
		--regions "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"

		"""

}


process FILTER_SNVS {
    publishDir "${params.outfolder}/${params.runID}/SNV/", mode: 'copy', overwrite: true
    tag "${sample}"
    label 'gatk'
	label 'mem_16GB'
	label 'core_8'
	input:
		tuple val(sample), path(vcf)
		path(fasta)

	output:
		tuple path("${sample}_deepvariant_filtered.vcf.gz"), path("${sample}_deepvariant_filtered.vcf.gz.tbi")
	script:
		"""

		bcftools view -f PASS ${vcf} -Ou | \
		bcftools norm -a -m - -f ${fasta} -Ou | \
		bcftools annotate --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' -Ou | \
		bcftools sort -Oz -o ${sample}_deepvariant_filtered.vcf.gz

		tabix -p vcf ${sample}_deepvariant_filtered.vcf.gz

		"""

}

process FILTER_AND_MERGE_SNVS {
    publishDir "${params.outfolder}/${params.runID}/SNV/", mode: 'copy', overwrite: false
    label 'glnexus'
	label 'mem_256GB'
	label 'core_36'
	input:
		path(vcf)
		path(fasta)
	output:
		tuple path("multisample_deepvariant_filtered.vcf.gz"), path("multisample_deepvariant_filtered.vcf.gz.tbi")
	script:
		"""

		glnexus_cli --config DeepVariant ${vcf} | \
		bcftools sort -Ou -m ${task.memory} | \
		bcftools annotate --set-id +'%CHROM\\_%POS\\_%REF\\_%ALT' -Ou | \
		bcftools view --threads ${task.cpus} -f PASS -Oz -o multisample_deepvariant_filtered.vcf.gz

		tabix -p vcf multisample_deepvariant_filtered.vcf.gz

		"""

}