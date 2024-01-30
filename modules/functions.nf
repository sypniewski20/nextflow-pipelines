def Read_samplesheet(samplesheet) {
	Channel
	    .fromPath (samplesheet)
	    .ifEmpty ( "Samplesheet empty." )
	    .splitCsv ( header:true, sep:'\t' )
	    .map { row -> [ row.sampleID, file(row.R1, checkIfExists: true), file(row.R2, checkIfExists: true) ]}
	    .set { ch_samples_initial }
}

def Read_bam_checkpoint(bam_checkpoint_sheet) {
	Channel
		.fromPath (bam_checkpoint_sheet)
	    .ifEmpty ( "Samplesheet empty." )
	    .splitCsv ( header:true, sep:'\t' )
	    .map { row -> [ row.sampleID, file(row.bam, checkIfExists: true), file(row.bai, checkIfExists: true) ]}
	    .set { ch_samples_checkpoint }
}

process FILTER_CNV_VCF {
	publishDir "${params.outfolder}/${params.runID}/CNV", pattern: "${sample}", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_8GB'
	label 'core_4'
	input:
		tuple val(sample), val(caller), path(vcf)
	output:
		tuple val(sample), path("${sample}_${caller}_cnv_sorted.vcf.gz"), path("${sample}_${caller}_cnv_sorted.vcf.gz.tbi")
	script:
		"""
            
		bcftools sort --threads ${task.cpus} ${delly} -Ov | \
		bcftools view -f "PASS" -Oz -o ${sample}_${caller}_cnv_sorted.vcf.gz

		tabix -pf vcf ${sample}_${caller}_cnv_sorted.vcf.gz
		"""
}
