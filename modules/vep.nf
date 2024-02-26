process FILTER_AND_MERGE_SNVS {
    publishDir "${params.outfolder}/${params.runID}/SNV/", mode: 'copy', overwrite: true
    label 'glnexus'
	label 'mem_256GB'
	label 'core_36'
	input:
		path(vcf)
		path(fasta)
	output:
		tuple path("multisample_deepvariant_filtered.vcf.gz"), path("multisample_deepvariant_filtered.vcf.gz.tbi")
	script:
		"""

        vep \
        --cache \
        --dir_cache /scratch/references/vep_cache \
        --species homo_sapiens \
        --fasta /scratch/references/fasta/GRCh38/GRCh38.fa \
        --assembly GRCh38 \
        --offline \
        --no_stats \
        --buffer_size 10000 \
        --compress_output bgzip \
        -i ${vcf} \
        -o $2 \
        --vcf \
        --fork ${task.cpus} \
        --force_overwrite \
        --pick \
        -e \
        --max_sv_size 100000000000000 \
        --custom file=${cnv_resource},short_name=gnomad,fields=PC%EVIDENCE%SVTYPE,format=vcf,type=within,reciprocal=1,overlap_cutoff=80 \
        --custom file=${clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN

		"""

}
