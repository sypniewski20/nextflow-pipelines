process MANTA_GERMLINE {
	publishDir "${params.outfolder}/${params.runID}/SV/manta", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'manta'
	label 'mem_8GB'
	label 'core_18'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
	output:
		tuple val(sample), path("manta/results/variants/diploidSV.vcf.gz")
	script:
		"""
            
        /manta/bin/configManta.py \
		--bam ${bam} \
		--referenceFasta ${fasta} \
		--runDir manta 

		manta/runWorkflow.py -j ${task.cpus}

		"""
}

process MANTA_EXOME_GERMLINE {
	publishDir "${params.outfolder}/${params.runID}/SV/manta", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'manta'
	label 'mem_8GB'
	label 'core_18'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
	output:
		tuple val(sample), path("manta/results/variants/diploidSV.vcf.gz")
	script:
		"""
            
        /manta/bin/configManta.py \
		--bam ${bam} \
		--referenceFasta ${fasta} \
		--runDir manta \
		--exome 

		./manta/runWorkflow.py -j ${task.cpus}

		"""
}

process MANTA_WES_JOINT {
	publishDir "${params.outfolder}/${params.runID}/SV", mode: 'copy', overwrite: true
	label 'manta'
	label 'mem_8GB'
	label 'core_18'
	input:
		path(bam)
		path(bai)
		path(fasta)
		path(fasta_fai)
	output:
		path("manta/results/variants/diploidSV.vcf.gz")
	script:
		"""

		BAM=\$( for i in \$( ls *.bam ); do echo "--bam \$i"; done )

        /manta/bin/configManta.py \
		--referenceFasta ${fasta} \
		--runDir manta \
		--exome \
		\$( echo \$BAM )

		./manta/runWorkflow.py -j ${task.cpus}

		"""
}

process MANTA_WGS_JOINT {
	publishDir "${params.outfolder}/${params.runID}/SV/manta", mode: 'copy', overwrite: true
	label 'manta'
	label 'mem_8GB'
	label 'core_18'
	input:
		path(bam)
		path(bai)
		path(fasta)
		path(fasta_fai)
	output:
		path("manta/results/variants/diploidSV.vcf.gz")
	script:
		"""
            
		BAM=\$( for i in \$( ls *.bam ); do echo "--bam \$i"; done )

        /manta/bin/configManta.py \
		--referenceFasta ${fasta} \
		--runDir manta \
		\$BAM

		./manta/runWorkflow.py -j ${task.cpus}

		"""
}

process MANTA_FILTER_VCF {
	publishDir "${params.outfolder}/${params.runID}/SV/manta", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_8GB'
	label 'core_4'
	input:
		tuple val(sample), path(manta)
	output:
		tuple path("${sample}_manta_sorted.vcf.gz"), path("${sample}_manta_sorted.vcf.gz.tbi")
	script:
		"""
            
		bcftools view -f PASS --threads ${task.cpus} ${manta} -Ou | \
		bcftools sort -Oz -o ${sample}_manta_sorted.vcf.gz

		tabix -p vcf ${sample}_manta_sorted.vcf.gz

		"""
}