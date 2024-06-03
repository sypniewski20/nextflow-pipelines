process DELLY_SV_CALL {
	tag "${sample}"
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
		path(contigs_bed)
	output:
		tuple val(sample), path("${sample}_delly_sv.vcf")
	script:
		"""

        delly call \
			  -g ${fasta} \
     		  -x ${contigs_bed} \
			  ${bam} > ${sample}_delly_sv.vcf

		"""
}

process DELLY_CNV_CALL {
	tag "${sample}"
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		tuple val(sample), path(delly_sv)
		path(fasta)
		path(fasta_fai)
		path(delly_map)
	output:
		tuple val(sample), path("${sample}_delly_cnv.vcf")
	script:
		"""

		delly cnv \
			  -o ${sample}_delly_cnv_temp.vcf \
			  -g ${fasta} \
			  -m ${delly_map} \
			  -l ${delly_sv} \
			  ${bam}


		delly classify -f germline -o ${sample}_delly_cnv.vcf ${sample}_delly_cnv_temp.vcf

		"""
}

process FILTER_DELLY {
    publishDir "${params.outfolder}/${params.runID}/SV/delly", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(vcf)
	output:
		tuple val(sample), path("${sample}_delly_sv_filtered.vcf.gz"), path("${sample}_delly_sv_filtered.vcf.gz.tbi")
	script:
		"""

		bcftools view -f PASS --threads ${task.cpus} ${vcf} -Ou | \
		bcftools sort -Oz -o ${sample}_delly_sv_filtered.vcf.gz

		tabix -p vcf ${sample}_delly_sv_filtered.vcf.gz

		"""
}


process DELLY_MERGE_SITES {
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		path(vcf)
	output:
		path("sites.vcf")
	script:
		"""

        delly merge ${vcf} -o sites.vcf

		"""
}

process DELLY_MERGED_SITES_CALL {
	tag "${sample}"
	label 'delly'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		file(fasta)
		path(contigs_bed)
		path(sites)
	output:
		tuple val(sample), path("${sample}_delly_joint.vcf")
	script:
		"""

        delly call \
			  -v ${sites} \
			  -x ${contigs_bed} \
			  -g ${fasta}/${fasta}.fa \
			  -o ${sample}_delly_joint.vcf
			  ${bam} 

		"""
}