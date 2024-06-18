process GRIDSS {
	publishDir "${params.outfolder}/${params.runID}/SV/gridss", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gridss'
	label 'mem_32GB'
	label 'core_8'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
    	path(fasta_fai)
        path(fasta_dict)
    	path(fasta_sa)
    	path(fasta_bwt)
    	path(fasta_ann)
    	path(fasta_amb)
    	path(fasta_pac)
	output:
		tuple val(sample), path("${sample}_gridss.vcf")
	script:
		"""

		gridss \
        -r ${fasta} \
		-o ${sample}_gridss.vcf \
		--skipsoftcliprealignment \
		-b ${params.exclude_bed} \
		${bam}

		"""
}

process GRIDSS_JOINT {
	label 'gridss'
	label 'mem_10GB'
	label 'core_4'
	input:
		path(bam), 
		path(bai),
		path(fasta)
    	path(fasta_fai)
        path(fasta_dict)
    	path(fasta_sa)
    	path(fasta_bwt)
    	path(fasta_ann)
    	path(fasta_amb)
    	path(fasta_pac)
	output:
		path("mutlisample_gridss.vcf")
	script:
		"""

		gridss \
        -r ${fasta} \
		--skipsoftcliprealignment \
		-o mutlisample_gridss.vcf \
		-b ${params.exclude_bed} \
		${bam}

		"""
}

process VIRUS_BREAKEND {
	publishDir "${params.outfolder}/${params.runID}/SV/gridss", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gridss'
	label 'mem_64GB'
	label 'core_4'
	input:
		tuple val(sample), path(bam)
	output:
		tuple val(sample), path("${sample}_virusbreakend.vcf")
	script:
		"""

        virusbreakend \
        -t ${task.cpus} \
        -o ${sample}_virusbreakend.vcf \
        --db ${params.virusbreakend_db} \
        ${bam}

		"""
}

process GRIDSS_FILTER_VCF {
	publishDir "${params.outfolder}/${params.runID}/SV/gridss", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_4'
	input:
		tuple val(sample), path(vcf)
	output:
		tuple val(sample), path("${vcf}_pass.gz"), path("${vcf}_pass.gz.tbi")
	script:
		"""

		bcftools view -f 'PASS' ${vcf} -Oz -o ${vcf}_pass.gz
        tabix -p ${vcf}_pass.gz

		"""
}