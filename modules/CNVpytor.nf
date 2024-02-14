process CNVPYTOR_CALL {
	tag "${sample}"
	label 'cnvpytor'
	label 'mem_16GB'
	label 'core_16'
	input:
		tuple val(sample), path(bam), path(bai)
	output:
		tuple val(sample), path("${sample}.vcf"), emit: vcf
		tuple val(sample), path("${sample}.pytor"), emit: pytor
	script:
		"""

		#!/usr/bin/python3

		import cnvpytor,os
		app = cnvpytor.Root("${sample}.pytor", create=True, max_cores=${task.cpus})
		app.rd(["${bam}"])
		app.calculate_histograms([${params.bins}])
		app.partition([${params.bins}])
		calls=app.call([${params.bins}])

		view = cnvpytor.Viewer(["${sample}.pytor"], params={} )
		view.bin_size = int(${params.bins})
		view.print_filename = "${sample}.vcf"
		view.print_calls_file()

		"""
}

process CNVPYTOR_FILTER_VCF {
	publishDir "${params.outfolder}/${params.runID}/CNV", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'	
	label 'mem_8GB'
	label 'core_4'
	input:
		tuple val(sample), path(vcf)
	output:
		tuple val(sample), path("${sample}_cnvpytor_sorted.vcf.gz"), path("${sample}_cnvpytor_sorted.vcf.gz.tbi")
	script:
		"""
            
		bcftools sort ${vcf} -Ov | \
		bcftools view -f PASS --threads ${task.cpus} -Oz -o ${sample}_cnvpytor_sorted.vcf.gz

		tabix -p vcf ${sample}_cnvpytor_sorted.vcf.gz
		"""
}


process CNVPYTOR_MERGE_CALLS {
	tag "${sample}"
	label 'cnvpytor'
	label 'mem_16GB'
	label 'core_36'
	input:
		path(pytor)
	output:
		path("multisample.vcf")
	script:
		"""

		#!/usr/bin/python3

		import cnvpytor,os

		view = cnvpytor.Viewer(["${pytor}"], params={} )
		view.Q0_range = [0, 0.5]
		view.size_range = int([100000, 'inf'])
		view.print_filename = "multisample.vcf"
		view.print_calls_file()

		"""
}

process CNVPYTOR_FILTER_MERGED {
	publishDir "${params.outfolder}/${params.runID}/CNV", mode: 'copy', overwrite: true
	label 'gatk'
	label 'mem_16GB'
	label 'core_36'
	input:
		path(vcf)
	output:
		tuple path("multisample_cnvpytor_sorted.vcf.gz"), path("multisample_cnvpytor_sorted.vcf.gz.tbi")
	script:
		"""
		
		bcftools sort ${vcf} -Ou | \
		bcftools view --threads ${task.cpus} -f PASS -Oz -o multisample_cnvpytor_sorted.vcf.gz

		tabix -p vcf multisample_cnvpytor_sorted.vcf.gz

		"""
}