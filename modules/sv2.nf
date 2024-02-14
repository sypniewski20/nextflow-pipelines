process SV2_REGENOTYPE {
	tag "${sample}"
	label 'sv2'	
	label 'mem_8GB'
	label 'core_4'
	input:
		tuple val(sample), path(bam), path(bai)
		tuple val(sample), path(delly), path(delly_tbi)
		tuple val(sample), path(manta), path(manta_tbi)
		tuple val(sample), path(cnvpytor), path(cnvpytor_tbi)
		tuple val(sample), path(erds), path(erds_tbi)
		tuple val(sample), path(snv), path(snv_tbi)
		path(ped_file)
		path(fasta)
	output:
        tuple val(sample), path("sv2_genotypes/${sample}_sv2.vcf")
	script:
		"""

		sv2 -hg38 ${fasta}/${fasta}.fa

        sv2 -g hg38 \
			-i ${bam} \
			-v ${delly} ${manta} ${cnvpytor} ${erds} \
			-snv ${snv} \
			-p ${ped_file} \
			-o ${sample}_sv2 \
			-L log.txt

		"""
}

process SV2_FILTER_VCF {
	publishDir "${params.outfolder}/${params.runID}/CNV", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_8GB'
	label 'core_4'
	input:
		tuple val(sample), path(sv2)
	output:
		tuple val(sample), path("${sample}_sv2_sorted.vcf.gz"), path("${sample}_sv2_sorted.vcf.gz.tbi"), emit: vcf
		tuple val(sample), path("${sample}_sv2.tsv"), emit: tsv
	script:
		"""

		bgzip -@ ${task.cpus} ${sv2}
		tabix -p vcf ${sv2}.gz

		bcftools sort ${sv2}.gz -Ov | \
		bcftools view -f PASS --threads ${task.cpus} -Oz -o ${sample}_sv2_sorted.vcf.gz

		tabix -p vcf ${sample}_sv2_sorted.vcf.gz

		echo -e "SAMPLE\tCHROM\tPOS\tID\tSVTYPE\tSVLEN\tCYTOBAND\tDESCRIPTION\tGT\tCN" > ${sample}_sv2.tsv

		bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%ID\t%SVTYPE\t%SVLEN\t%CYTOBAND\t%DESCRIPTION\t%GT\t%CN\n]' ${sample}_sv2_sorted.vcf.gz >> ${sample}_sv2.tsv

		"""
}
