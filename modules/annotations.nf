process ANNOT_SV {
    publishDir "${params.outfolder}/${params.runID}/SV", mode: 'copy', overwrite: true
	label 'annotsv'
	input:
		tuple file(vcf), file(tbi)
	output:
		path("annotated_${vcf}.tsv")
	script:
		"""

		AnnotSV -SVinputFile ${vcf} -genomeBuild ${params.genome} -outputFile annotated_${vcf} -outputDir .

		"""
}	

process VEP {
    publishDir "${params.outfolder}/${params.runID}/SNV/", mode: 'copy', overwrite: true
    label 'vep'
	label 'mem_64GB'
	label 'core_36'
	input:
		path(vcf)
		path(fasta)
	output:
		path("${vcf}.vep.tsv.gz")
	script:
		"""

        vep \
        --cache \
        --dir_cache /scratch/references/vep_cache \
        --species homo_sapiens \
        --fasta ${fasta}/${fasta}.fa \
        --assembly ${params.genome} \
        --offline \
        --no_stats \
        --buffer_size 10000 \
        --compress_output bgzip \
        -i ${vcf} \
        -o ${vcf}.vep.tsv.gz \
        --tab \
        --fork ${task.cpus} \
        --force_overwrite \
        --pick \
        -e \
        --max_sv_size 100000000000000 \
        --custom file=${params.clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN

		"""

}
