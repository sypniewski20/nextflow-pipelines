process DELLY_CNV_CALL {
	tag "${sample}"
	label 'delly'
	label 'mem_16GB'
	label 'core_16'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(centromeres)
	output:
		tuple val(sample), path("${sample}_delly.vcf.gz")
	script:
		"""

		export OMP_NUM_THREADS=${task.cpus}

        delly call \
			  -x ${centromeres} \
			  -g ${fasta}/${fasta}.fa \
			  ${bam} | \
		bgzip -c -@ ${task.cpus} > ${sample}_delly.vcf.gz

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

process DELLY_MERGE_SITES {
	label 'delly'
	label 'mem_16GB'
	label 'core_16'
	input:
		tuple val(sample), path(vcf)
	output:
		path("sites.vcf.gz")
	script:
		"""

		export OMP_NUM_THREADS=${task.cpus}

        delly merge ${vcf} | \
		bgzip -c -@ ${task.cpus} > sites.vcf.gz

		"""
}

process DELLY_MERGED_SITES_CALL {
	tag "${sample}"
	label 'delly'
	label 'mem_16GB'
	label 'core_16'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(centromeres)
		path(sites)
	output:
		tuple val(sample), path("${sample}_delly.vcf.gz")
	script:
		"""

		export OMP_NUM_THREADS=${task.cpus}

        delly call \
			  -v ${sites} \
			  -x ${centromeres} \
			  -g ${fasta}/${fasta}.fa \
			  ${bam} | \
		bgzip -c -@ ${task.cpus} > ${sample}_delly.vcf.gz

		"""
}

process DELLY_MERGED_FILTER {
	label 'delly'
	label 'mem_16GB'
	label 'core_16'
	input:
		path(vcf)
	output:
		tuple path("multisample_delly_filtered.vcf.gz"), path("multisample_delly_filtered.vcf.gz")
	script:
		"""

		export OMP_NUM_THREADS=${task.cpus}

		bcftools merge ${vcf} -Ou | \
        delly filter \
			  -f germline | \
		bgzip -c -@ ${task.cpus} > multisample_delly_filtered.vcf.gz

		tabix -p multisample_delly_filtered.vcf.gz

		"""
}