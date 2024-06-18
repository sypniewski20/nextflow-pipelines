process DELLY_SV_CALL {
	publishDir "${params.outfolder}/${params.runID}/SV/delly", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
	output:
		tuple val(sample), path("${sample}_delly_sv.vcf.gz")
	script:
		"""

        delly call \
			  -g ${fasta} \
     		  -x ${params.exclude_bed} \
			  ${bam} | bgzip -c > ${sample}_delly_sv.vcf.gz

		"""
}

process DELLY_MERGE_SV_SITES {
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		path(vcf)
	output:
		path("sites.bcf")
	script:
		"""

        delly merge ${vcf} -o sites.bcf

		"""
}

process DELLY_JOINT_SV_CALL {
	tag "${sample}"
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
		path(sites)
	output:
		tuple val(sample), path("${sample}_delly_joint.bcf"), path("${sample}_delly_joint.bcf.csi")
	script:
		"""

        delly call \
			  -v ${sites} \
			  -x ${params.exclude_bed} \
			  -g ${fasta} \
			  -o ${sample}_delly_joint.bcf \
			  ${bam} 

		"""
}

process DELLY_SV_MERGE {
    publishDir "${params.outfolder}/${params.runID}/SV/delly", mode: 'copy', overwrite: true
	label 'gatk'
	label 'mem_16GB'
	label 'core_1'
	input:
		path(bcf)
		path(csi)
	output:
		path("delly_sv_merged.vcf.gz"), emit: vcf
		path("delly_sv_merged.vcf.gz"), emit: tbi
	script:
		"""

		bcftools merge ${bcf} -Oz -o delly_sv_merged.vcf.gz
		tabix -p vcf delly_sv_merged.vcf.gz

		"""
}

process DELLY_CNV_CALL {
	publishDir "${params.outfolder}/${params.runID}/SV/delly", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
		file(sites)
	output:
		tuple val(sample), path("${sample}.bcf"), path("${sample}.bcf.csi")
	script:
		"""

		delly cnv -i 100000 \
				  -j 100000 \
				  -w 100000 \
				  -o ${sample}.bcf \
				  -g ${fasta} \
				  -m ${params.delly_map} \
				  -l ${sites} \
				  ${bam}

		"""
}

process DELLY_MERGE_CNV_SITES {
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		path(bcf)
		path(csi)
	output:
		path("sites_cnv.bcf")
	script:
		"""

		delly merge -e \
					-p \
					-o sites_cnv.bcf \
					-m 1000 \
					-n 100000 \
					${bcf}

		"""
}

process DELLY_JOINT_CNV_CALL {
	tag "${sample}"
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
		path(sites)
	output:
		tuple val(sample), path("${sample}_delly_cnv.bcf"), path("${sample}_delly_cnv.bcf.csi")
	script:
		"""

		delly cnv -i 100000 \
				  -j 100000 \
				  -w 100000 \
				  -u \
				  -v ${sites} \
				  -g ${fasta} \
				  -m ${params.delly_map} \
				  -o  ${sample}_delly_cnv.bcf \
				  ${bam}

		"""
}

process DELLY_CNV_MERGE {
    publishDir "${params.outfolder}/${params.runID}/SV/delly", mode: 'copy', overwrite: true
	label 'gatk'
	label 'mem_16GB'
	label 'core_1'
	input:
		path(bcf)
		path(csi)
	output:
		path("delly_cnv_merged.vcf.gz"), emit: vcf
		path("delly_cnv_merged.vcf.gz.tbi"), emit: tbi
	script:
		"""

		bcftools merge ${bcf} -Oz -o delly_cnv_merged.vcf.gz
		tabix -p vcf delly_cnv_merged.vcf.gz

		"""
}

process DELLY_MERGE_SV_CNV_CALLS {
    publishDir "${params.outfolder}/${params.runID}/SV/delly", mode: 'copy', overwrite: true
	label 'gatk'
	label 'mem_16GB'
	label 'core_1'
	input:
		path(sv)
		path(sv_tbi)
		path(cnv)
		path(cnv_tbi)
	output:
		tuple path("multisample_delly_calls.vcf.gz"), path("multisample_delly_calls.vcf.gz.tbi")
	script:
		"""

		bcftools concat ${sv} ${cnv} -Oz -o multisample_delly_calls.vcf.gz
		tabix -p vcf multisample_delly_calls.vcf.gz

		"""
}