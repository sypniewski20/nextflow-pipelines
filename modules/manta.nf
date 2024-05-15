process MANTA_GERMLINE {
	tag "${sample}"
	label 'manta'
	label 'mem_8GB'
	label 'core_18'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
		path(contigs)
		path(contigs_tbi)
	output:
		tuple val(sample), path("manta/results/variants/diploidSV.vcf.gz")
	script:
		"""
            
        /manta/bin/configManta.py \
		--bam ${bam} \
		--referenceFasta ${fasta} \
		--runDir manta \
		--callRegions ${contigs}

		manta/runWorkflow.py -j ${task.cpus}

		"""
}

process MANTA_EXOME_GERMLINE {
	tag "${sample}"
	label 'manta'
	label 'mem_8GB'
	label 'core_18'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
		path(contigs)
		path(contigs_tbi)
	output:
		tuple val(sample), path("manta/results/variants/diploidSV.vcf.gz")
	script:
		"""
            
        /manta/bin/configManta.py \
		--bam ${bam} \
		--referenceFasta ${fasta} \
		--runDir manta \
		--exome \
		--callRegions ${contigs}

		manta/runWorkflow.py -j ${task.cpus}

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
		tuple path("${sample}_manta_cnv_sorted.vcf.gz"), path("${sample}_manta_cnv_sorted.vcf.gz.tbi")
	script:
		"""
            
		bcftools view -f PASS --threads ${task.cpus} ${manta} -Ou | \
		bcftools sort -Oz -o ${sample}_manta_cnv_sorted.vcf.gz

		tabix -p vcf ${sample}_manta_cnv_sorted.vcf.gz

		"""
}

process MANTA_MERGE_VCF {
	publishDir "${params.outfolder}/${params.runID}/SV/manta", mode: 'copy', overwrite: true
	label 'gatk'
	label 'mem_8GB'
	label 'core_4'
	input:
		path(vcf)
		path(tbi)
	output:
		tuple path("multisample_manta_cnv_sorted.vcf.gz"), path("multisample_manta_cnv_sorted.vcf.gz.tbi")
	script:
		"""
            
		bcftools merge --missing-to-ref ${vcf} -Ov | \
		bcftools sort -Oz -o multisample_manta_cnv_sorted.vcf.gz

		tabix -p vcf multisample_manta_cnv_sorted.vcf.gz
		
		"""
}