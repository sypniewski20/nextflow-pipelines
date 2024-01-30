process DOWNLOAD_FASTA {
    publishDir "${params.references}/fasta/${fasta_name}", mode: 'copy', overwrite: true
    label 'wgs_tools'
	label 'mem_16GB'
	label 'core_16'
	output:
		tuple val("${fasta_name}"), path("${fasta_name}.fa")
	script:
		"""

        ${fasta_name}=\$( basename ${params.fasta_url} | sed 's/\\.fa\\.gz//g' )

        wget ${params.fasta_url}
        bgzip -@ ${task.cpus} -d \$( echo ${fasta_name}.fa )

		"""
}

process FAI_FASTA {
    publishDir "${params.references}/fasta/${fasta_name}", mode: 'copy', overwrite: true
    label 'wgs_tools'
	label 'mem_8GB'
	label 'core_1'
    input:
        tuple val(fasta_name), file(fasta)
	output:
		tuple val("${fasta_name}"), file("${fasta}"), file("${fasta}.fai")
	script:
		"""

        samtools faidx ${fasta}

		"""
}

process BWA_INDEX {
    publishDir "${params.references}/fasta/${fasta_name}", mode: 'copy', overwrite: true
    label 'wgs_tools'
	label 'mem_48GB'
	label 'core_1'
    input:
        tuple val(fasta_name), file(fasta), file(fai)
	output:
		tuple val("${fasta_name}"), file("${fasta}"), file("${fasta}.fai")
	script:
		"""

        bwa index ${fasta}

		"""
}

process SEQ_DICT {
    publishDir "${params.references}/fasta/${fasta_name}", mode: 'copy', overwrite: true
    label 'wgs_tools'
	label 'mem_32GB'
	label 'core_1'
    input:
        tuple val(fasta_name), file(fasta), file(fai)
	output:
		file("${fasta}.dict")
	script:
		"""

         gatk CreateSequenceDictionary \
        -R ${fasta}

		"""
}

process BWA_IMAGE {
    publishDir "${params.references}/fasta/${fasta_name}", mode: 'copy', overwrite: true
    label 'wgs_tools'
	label 'mem_48GB'
	label 'core_1'
    input:
        tuple val(fasta_name), file(fasta), file(fai)
	output:
		file("${fasta}.img")
	script:
		"""

         gatk BwaMemIndexImageCreator \
        -I ${fasta} \
        -O ${fasta}.img

		"""
}