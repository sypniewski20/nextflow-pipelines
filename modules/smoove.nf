process SMOOVE {
	publishDir "${params.outfolder}/${params.runID}/SV/smoove", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'smoove'
	label 'mem_16GB'
	label 'core_8'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(fasta_fai)
	output:
		tuple val(sample), path("${sample}-smoove.genotyped.vcf.gz"), path("${sample}-smoove.genotyped.vcf.gz.csi")
	script:
		"""
            
        smoove call -x \
					-duphold \
					--genotype \
                    --name ${sample} \
                    --exclude /data/references/bed/GRCh38/exclude.cnvnator_100bp.GRCh38.20170403.bed \
                    --fasta ${fasta} \
                    -p ${task.cpus} \
					--outdir . \
					${bam}

		"""
}

process SMOOVE_ANNOTATE {
	publishDir "${params.outfolder}/${params.runID}/SV/smoove", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'smoove'
	label 'mem_8GB'
	label 'core_1'
	input:
		tuple val(sample), path(vcf), path(tbi)
	output:
		tuple val(sample), path("anno_${vcf}"), path("anno_${vcf}.tbi")
	script:
		"""
            
		smoove annotate --gff ${params.gff} ${vcf} | bgzip -c > anno_${vcf}
		tabix -p vcf anno_${vcf}

		"""
}



process SMOOVE_MERGE {
	publishDir "${params.outfolder}/${params.runID}/SV/smoove", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'smoove'
	label 'mem_8GB'
	label 'core_1'
	input:
		path(vcf)
		path(tbi)
		path(fasta)
	output:
		tuple path("${params.runID}_merged.sites.vcf.gz"), path("${params.runID}_merged.sites.vcf.gz.tbi")
	script:
		"""
            
        smoove merge --name ${params.runID}_merged \
					 --fasta ${fasta} \
					 --outdir . \
					 ${vcf}

		"""
}

process SMOOVE_JOINT_GENOTYPE {
	publishDir "${params.outfolder}/${params.runID}/SV/smoove", mode: 'copy', overwrite: true
	label 'smoove'
	label 'mem_16GB'
	label 'core_1'
	input:
		path(bam)
		path(bai)
		path(fasta)
		path(exclude_bed)
	output:
		tuple path("multisample-smoove.genotyped.vcf.gz"), path("multisample-smoove.genotyped.vcf.gz.tbi")
	script:
		"""
            
		smoove call -x \
					--name multisample \
					--exclude ${exclude_bed} \
					--fasta ${fasta} \
					-p ${task.cpus} \
					--genotype ${bam}

		tabix -p vcf multisample-smoove.genotyped.vcf.gz


		"""
}