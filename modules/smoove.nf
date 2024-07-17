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
                    --exclude ${params.exclude_bed} \
                    --fasta ${fasta} \
                    -p ${task.cpus} \
					--outdir . \
					${bam}

		"""
}

process SMOOVE_JOINT {
	publishDir "${params.outfolder}/${params.runID}/SV/smoove", mode: 'copy', overwrite: true
	label 'smoove'
	label 'mem_16GB'
	label 'core_8'
	input:
		path(bam), 
		path(bai),
		path(fasta)
		path(fasta_fai)
	output:
		tuple path("multisample-smoove.genotyped.vcf.gz"), path("multisample-smoove.genotyped.vcf.gz.csi")
	script:
		"""
            
        smoove call -x \
					-duphold \
					--genotype \
                    --name ${sample} \
                    --exclude ${params.exclude_bed} \
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
		path(fasta_fai)
	output:
		tuple val(sample), path("multisample-smoove.genotyped.vcf.gz"), emit: vcf
		tuple val(sample), path("multisample-smoove.genotyped.vcf.gz.tbi"), emit: tbi
	script:
		"""
            
		smoove call -x \
					--name multisample \
					--exclude ${params.exclude_bed} \
					--fasta ${fasta} \
					-p ${task.cpus} \
					--genotype ${bam}

		tabix -p vcf multisample-smoove.genotyped.vcf.gz


		"""
}


process SMOOVE_FILTER_VCF {
	publishDir "${params.outfolder}/${params.runID}/SV/smoove", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_8GB'
	label 'core_4'
	input:
		tuple val(sample), path(vcf), path(tbi)
	output:
		tuple val(sample), path("${sample}_smoove_sorted.vcf.gz"), path("${sample}_smoove_sorted.vcf.gz.tbi")
	script:
		"""
            
		bcftools view -f PASS --threads ${task.cpus} ${vcf} -Ou | \
		bcftools sort -Oz -o ${sample}_smoove_sorted.vcf.gz

		tabix -p vcf ${sample}_smoove_sorted.vcf.gz

		"""
}