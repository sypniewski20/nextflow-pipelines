process DELLY_SV_CALL {
	tag "${sample}"
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
	output:
		tuple val(sample), path("${sample}_delly_sv.bcf")
	script:
		"""

        delly call \
			  -g ${fasta} \
     		  -x ${params.exclude_bed} \
			  ${bam} > ${sample}_delly_sv.bcf

		"""
}

process DELLY_MERGE_SV_SITES {
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		path(vcf)
	output:
		path("sites.vcf")
	script:
		"""

        delly merge ${vcf} -o sites.bcf

		"""
}

process DELLY_MERGED_SV_SITES_CALL {
	tag "${sample}"
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		file(fasta)
		file(fasta_fai)
		path(sites)
	output:
		tuple val(sample), path("${sample}_delly_joint.bcf")
	script:
		"""

        delly call \
			  -v ${sites} \
			  -x ${params.exclude_bed} \
			  -g ${fasta} \
			  -o ${sample}_delly_joint.bcf
			  ${bam} 

		"""
}

process DELLY_SV_MERGE {
    publishDir "${params.outfolder}/${params.runID}/SV/delly", mode: 'copy', overwrite: true
	label 'gatk'
	label 'mem_16GB'
	label 'core_1'
	input:
		path(vcf)
	output:
		tuple path("delly_sv_merged.vcf.gz"), path("delly_sv_merged.vcf.gz.tbi")
	script:
		"""

		bcftools merge ${vcf} -Oz -o delly_sv_merged.vcf.gz
		tabix -p vcf delly_sv_merged.vcf.gz

		"""
}

process DELLY_CNV_CALL {
	tag "${sample}"
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		file(fasta)
		file(fasta_fai)
		tuple file(sites), file(sites_tbi)
	output:
		tuple val(sample), path("${sample}_delly_cnv.bcf")
	script:
		"""

		delly cnv -o ${sample}.bcf \
				  -g ${fasta} \
				  -m ${params.delly_map} \
				  -l ${sites} \
				  ${bam}

		"""
}

process DELLY_CNV_MERGE_SITES {
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		file(vcf)
	output:
		path("sites_cnv.bcf")
	script:
		"""

		delly merge -e \
					-p \
					-o ${sites} \
					-m 1000 \
					-n 100000 \
					${vcf}

		"""
}

process DELLY_JOINT_CNV_CALL {
	tag "${sample}"
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		file(fasta)
		file(fasta_fai)
		file(sites)
	output:
		tuple val(sample), path("${sample}_delly_cnv.bcf")
	script:
		"""

		delly cnv -u \
				  -v ${sites} \
				  -g ${fasta} \
				  -m ${params.delly_map} \
				  -o  ${sample}_delly_cnv.bcf \
				  input.bam

		"""
}

process DELLY_CNV_MERGE_SAMPLES {
    publishDir "${params.outfolder}/${params.runID}/SV/delly", mode: 'copy', overwrite: true
	label 'gatk'
	label 'mem_16GB'
	label 'core_1'
	input:
		file(vcf)
	output:
		tuple val(sample), path("${sample}_delly_cnv.bcf")
	script:
		"""

		bcftools merge ${vcf} -Oz -o delly_sv_merged.vcf.gz
		tabix -p vcf delly_cnv_merged.vcf.gz

		"""
}