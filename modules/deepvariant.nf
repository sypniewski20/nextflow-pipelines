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
		path(contigs_bed)
	output:
		tuple val(sample), path("${sample}_deepvariant.vcf.gz")
	script:
		"""

        /opt/deepvariant/bin/run_deepvariant \
        --model_type WGS \
        --ref ${fasta} \
        --reads ${bam} \
		--output_vcf ${sample}_deepvariant.vcf.gz \
        --num_shards ${task.cpus} \
		--regions ${contigs_bed}

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
		path(contigs_bed)
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
		--regions chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX

		"""

}

process DEEP_VARIANT_WGS_GVCF {
    tag "${sample}"
    label 'deepvariant'
	label 'mem_96GB'
	label 'core_32'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
		path(contigs_bed)
	output:
		tuple val(sample), path("${sample}_deepvariant.g.vcf.gz")
	script:
		"""

        /opt/deepvariant/bin/run_deepvariant \
        --model_type WGS \
        --ref ${fasta} \
        --reads ${bam} \
        --output_gvcf ${sample}_deepvariant.g.vcf.gz \
		--output_vcf ${sample}_deepvariant.vcf.gz \
        --num_shards ${task.cpus} \
		--regions ${contigs_bed}

		"""

}

process DEEP_VARIANT_WES_GVCF {
    tag "${sample}"
    label 'deepvariant'
	label 'mem_96GB'
	label 'core_32'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
		path(contigs_bed)
	output:
		tuple val(sample), path("${sample}_deepvariant.g.vcf.gz")
	script:
		"""

        /opt/deepvariant/bin/run_deepvariant \
        --model_type WES \
        --ref ${fasta} \
        --reads ${bam} \
        --output_gvcf ${sample}_deepvariant.g.vcf.gz \
		--output_vcf ${sample}_deepvariant.vcf.gz \
        --num_shards ${task.cpus} \
		--regions ${contigs_bed}

		"""

}

process FILTER_SNVS {
    publishDir "${params.outfolder}/${params.runID}/SNV/deep_variant", mode: 'copy', overwrite: true
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

process GLNEXUS_WES {
    label 'glnexus'
	label 'mem_64GB'
	label 'core_36'
	input:
		path(vcf)
	output:
		path("multisample_deepvariant.vcf.gz")
	script:
		"""

		glnexus_cli \
		--threads ${task.cpus} \
		--config DeepVariantWES \
		--bed "${params.contigs_bed}" \
		${vcf} | \
		bcftools view -Oz -o multisample_deepvariant.vcf.gz

		"""

}

process NORM_MULTISAMPLE {
    publishDir "${params.outfolder}/${params.runID}/SNV/deep_variant", mode: 'copy', overwrite: true
    label 'gatk'
	label 'mem_16GB'
	label 'core_16'
	input:
		path(vcf)
		path(fasta)
	output:
		tuple path("norm_${vcf}"), path("norm_${vcf}.tbi")
	script:
		"""

		bcftools +fill-tags ${vcf} -Ou -- -t AF,AC,AN | \
		bcftools norm -f ${fasta} -a -m - -Ou | \
		bcftools sort -Oz -o norm_${vcf}

		tabix -p vcf norm_${vcf}

		"""

}
