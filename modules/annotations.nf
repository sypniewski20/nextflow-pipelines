process ANNOT_SV {
    publishDir "${params.outfolder}/${params.runID}/SV", mode: 'copy', overwrite: true
	label 'annotsv'
	input:
		tuple file(vcf), file(tbi)
	output:
		path("annotated_${vcf}.tsv")
	script:
		"""

		AnnotSV -SVinputFile ${vcf} \
                -genomeBuild ${params.genome} \
                -outputFile annotated_${vcf} \
                -outputDir .

		"""
}	

process VEP {
    publishDir "${params.outfolder}/${params.runID}/SNV/", mode: 'copy', overwrite: true
    label 'vep'
	label 'mem_64GB'
	label 'core_36'
	input:
		tuple file(vcf), file(tbi)
		path(fasta)
        path(fasta_fai)
	output:
		path("${vcf}.vep.tsv.gz")
	script:
		"""

        vep \
        --cache \
        --dir_cache ${params.vep_cache} \
        --species homo_sapiens \
        --fasta ${fasta} \
        --assembly ${params.genome} \
        --buffer_size 1000000 \
        -i ${vcf} \
        -o ${vcf}.vep.tsv.gz \
        --format vcf \
        --compress_output bgzip \
        --tab \
        --fork ${task.cpus} \
        --force_overwrite \
        --e \
        --custom file=${params.clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \
        --plugin REVEL,file=${params.REVEL}


		"""

}

process VEP_SV {
    publishDir "${params.outfolder}/${params.runID}/SNV/", pattern: "${vcf}.vep.tsv.gz*", mode: 'copy', overwrite: true
    tag "${sample}"
    label 'vep'
	label 'mem_64GB'
	label 'core_36'
	input:
		tuple val(sample), file(vcf), file(tbi)
		path(fasta)
        path(fasta_fai)
	output:
		path("${vcf}.vep.tsv.gz")
	script:
		"""

        vep \
        --cache \
        --dir_cache ${params.vep_cache} \
        --species homo_sapiens \
        --fasta ${fasta} \
        --assembly ${params.genome} \
        --buffer_size 1000000 \
        --format vcf \
        --compress_output bgzip \
        -i ${vcf} \
        -o ${vcf}.vep.tsv.gz \
        --tab \
        --fork ${task.cpus} \
        --force_overwrite \
        --e \
        --plugin StructuralVariantOverlap,file=/data/references/germline_resource/gnomAD4/gnomad_sv_v4.0_complete.vcf.gz \
        --max_sv_size 1000000000


		"""

}


process SURVIVOR {
    publishDir "${params.outfolder}/${params.runID}/SV/SURVIVOR", mode: 'copy', overwrite: true
    tag "${sample}"
	label 'annotsv'
	input:
		tuple val(sample), file(manta), file(tbi)
        tuple val(sample), file(delly), file(tbi)
        tuple val(sample), file(smoove), file(tbi)
	output:
		tuple val(sample), path("${sample}_survivor_merged.vcf.gz"), path("${sample}_survivor_merged.vcf.gz.tbi")
	script:
		"""

        for vcf in *.vcf.gz; do 
         bgzip -d \$vcf 
        done      

        realpath *.vcf > sample_files

        SURVIVOR merge sample_files 1000 2 1 1 0 30 ${sample}_survivor_merged.vcf

        bcftools sort ${sample}_survivor_merged.vcf -Oz -o ${sample}_survivor_merged.vcf.gz
        tabix -p vcf ${sample}_survivor_merged.vcf.gz

		"""
}

process ANNOT_SV {
    publishDir "${params.outfolder}/${params.runID}/SV", mode: 'copy', overwrite: true
    tag "${sample}"
	label 'annotsv'
	input:
		tuple val(sample), file(vcf), file(tbi)
	output:
		path("annotated_${vcf}.tsv")
	script:
		"""

		AnnotSV -SVinputFile ${vcf} \
                -genomeBuild ${params.genome} \
                -outputFile annotated_${vcf} \
                -outputDir .

		"""
}