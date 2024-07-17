
process MUTECT2_NORMAL {
    publishDir "${params.outfolder}/${params.runID}/VCF/mutect2/germline", mode: 'copy', overwrite: true
	label 'gatk'
	label 'mem_32GB'
	label 'core_4'
	tag "${sampleID}"
	input:
		tuple val(sampleID), val(is_germline), file(bam), file(bai)
		path(fasta)
		path(exome_bed)
		val(interval_padding)
	output:
		path("${sampleID}_mutect2_germline.vcf.gz*")
	when:
	"${is_germline}" == "1"
	script:
		"""

		ln -sf \$( echo "\$( realpath ${exome_bed} ).tbi" ) .

		export OMP_NUM_THREADS=${task.cpus}

        gatk --java-options "-XX:ParallelGCThreads=${task.cpus}" Mutect2 \
		--max-mnp-distance 0 \
		-L ${exome_bed} \
		--interval-padding ${interval_padding} \
        -R ${fasta}/${fasta}.fa \
        -I ${bam} \
        -O ${sampleID}_mutect2_germline.vcf.gz \
		--native-pair-hmm-threads ${task.cpus}

		tabix -fp vcf ${sampleID}_mutect2_germline.vcf.gz
		"""
}

process MUTECT2_PON {
	label 'gatk'
	label 'mem_32GB'
	label 'core_8'
	input:
		file(vcf)
		path(fasta)
		path(exome_bed)
		val(interval_padding)
	output:
		path("pon_${params.runID}.vcf.gz")
	script:
		"""
		ls *vcf.gz > vcf.list

		ln -sf \$( echo "\$( realpath ${exome_bed} ).tbi" ) .

		gatk GenomicsDBImport -R ${fasta}/${fasta}.fa \
		--genomicsdb-workspace-path pon_db \
		-V vcf.list \
		-L ${exome_bed} \
		--interval-padding ${interval_padding} \
		--merge-input-intervals \
		--reader-threads ${task.cpus}

		gatk CreateSomaticPanelOfNormals -R ${fasta}/${fasta}.fa \
		-V gendb://pon_db \
		-O pon_${params.runID}.vcf.gz

		tabix -fp vcf pon_${params.runID}.vcf.gz

		"""
}

process MUTECT2_ROA_MODEL {
	label 'gatk'
	label 'mem_32GB'
	label 'core_4'
	tag "${sampleID}"
	input:
		tuple val(sampleID), val(is_germline), file(bam), file(bai)
        path(fasta)
        path(snv_resource)
		path(exome_bed)
		path(pon)
		val(interval_padding)
	output:
		tuple val(sampleID), val(is_germline), path("${sampleID}_unfiltered.vcf.gz"), path("${sampleID}_unfiltered.vcf.gz.tbi"), path("${sampleID}_unfiltered.vcf.gz.stats"), path("read-orientation-model.tar.gz"), emit: roa_model
	when:
	"${is_germline}" == "0"
	script:
		"""
		ln -s \$( echo "\$( realpath ${snv_resource} ).tbi" ) .
		ln -s \$( echo "\$( realpath ${pon} ).tbi" ) .
		ln -sf \$( echo "\$( realpath ${exome_bed} ).tbi" ) .

		export OMP_NUM_THREADS=${task.cpus}

        gatk --java-options "-XX:ParallelGCThreads=${task.cpus}" Mutect2 \
		--native-pair-hmm-threads ${task.cpus} \
        -R ${fasta}/${fasta}.fa \
        -I ${bam} \
		-L ${exome_bed} \
		--interval-padding ${interval_padding} \
        --germline-resource ${snv_resource} \
		--panel-of-normals ${pon} \
        --f1r2-tar-gz f1r2.tar.gz \
        -O ${sampleID}_unfiltered.vcf.gz

		gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz

		"""
}

process MUTECT2_PILEUP {
	label 'gatk'
	label 'mem_32GB'
	label 'core_4'
	tag "${sampleID}"
	input:
		tuple val(sampleID), val(is_germline), file(bam), file(bai)
        path(snv_resource)
		path(exome_bed)
		val(interval_padding)
	output:
		tuple val(sampleID), val(is_germline), path("getpileupsummaries.table")
	when:
	"${is_germline}" == "0"
	script:
		"""
		ln -s \$( echo "\$( realpath ${snv_resource} ).tbi" ) .
		ln -sf \$( echo "\$( realpath ${exome_bed} ).tbi" ) .


		gatk GetPileupSummaries \
		-I ${bam} \
		-L ${exome_bed} \
		--interval-padding ${interval_padding} \
		-V ${snv_resource} \
		-O getpileupsummaries.table

		"""
}

process MUTECT2_CONTAMINATION {
	label 'gatk'
	label 'mem_32GB'
	label 'core_4'
	tag "${sampleID}"
	input:
		tuple val(sampleID), val(is_germline), file(pileup)
	output:
		tuple val("${sampleID}"), path("segments.table"), path("calculatecontamination.table")
	when:
	"${is_germline}" == "0"
	script:
		"""
		gatk CalculateContamination \
		-I ${pileup} \
		-tumor-segmentation segments.table \
		-O calculatecontamination.table

		"""
}

def Fun_to_contamination_prep(roa_output, contamination) {
	roa_output
		.join( contamination )
}

process MUTECT2_ROA_FILTER {
	label 'gatk'
	label 'mem_32GB'
	label 'core_4'
	tag "${sampleID}"
	input:
		tuple val(sampleID), val(is_germline), file(vcf), file(tbi), file(stats), file(roa_model), file(segments_table), file(contamination_table)
		path(fasta)
	output:
		tuple val("${sampleID}"), path("${sampleID}_mutect2_ROA.vcf.gz"), path("${sampleID}_mutect2_ROA.vcf.gz.tbi")
	when:
	"${is_germline}" == "0"
	script:
		"""

		gatk FilterMutectCalls -V ${vcf} \
        --tumor-segmentation ${segments_table} \
		-R ${fasta}/${fasta}.fa \
        --contamination-table ${contamination_table} \
        --ob-priors ${roa_model} \
        -O ${sampleID}_mutect2_ROA.vcf.gz 

		tabix -fp vcf ${sampleID}_mutect2_ROA.vcf.gz

		"""
}

process MUTECT2_AS_FILTER_HEADER {
    publishDir "${params.outfolder}/${params.runID}/VCF/mutect2/somatic/", mode: 'copy', overwrite: true
	label 'gatk'
	label 'mem_32GB'
	label 'core_8'
	tag "${sampleID}"
	input:
		tuple val(sampleID), file(vcf), file(tbi)
		path(fasta)
	output:
		path("${sampleID}_mutect2_ROA_FAA_PASS.vcf.gz"), emit: vcf 
		path("${sampleID}_mutect2_ROA_FAA_PASS.vcf.gz.tbi"), emit: tbi
	script:
		"""

		bgzip -@ ${task.cpus} -d ${vcf}

		bcftools view -h ${sampleID}_mutect2_ROA.vcf | \
		perl -nle "if (/ID=AS_FilterStatus,/){ s/Number=A/Number=./ } print" > my.header.txt

		bcftools reheader -h my.header.txt -o ${sampleID}_mutect2_reheader_ROA.vcf ${sampleID}_mutect2_ROA.vcf 
		bgzip -@ ${task.cpus} ${sampleID}_mutect2_reheader_ROA.vcf

		tabix -fp vcf ${sampleID}_mutect2_reheader_ROA.vcf.gz

		bcftools view -f PASS ${sampleID}_mutect2_reheader_ROA.vcf.gz -Oz -o ${sampleID}_mutect2_ROA_FAA_PASS.vcf.gz

		tabix -fp vcf ${sampleID}_mutect2_ROA_FAA_PASS.vcf.gz

		"""
}

process MERGE_VCF {
	label 'gatk'
	label 'mem_32GB'
	label 'core_16'
    publishDir "${params.outfolder}/${params.runID}/SNV/mutect2/somatic/", pattern: "mutlisampleID_${params.runID}_mutect2_ROA_FAA_PASS.vcf.gz*", mode: 'copy', overwrite: true
	input:
		file(vcf)
		file(tbi)
	output:
		tuple path("mutlisampleID_${params.runID}_mutect2_ROA_FAA_PASS.vcf.gz"), path("mutlisampleID_${params.runID}_mutect2_ROA_FAA_PASS.vcf.gz.tbi"), emit: ch_somatic_output
	script:
		"""

		bcftools merge ${vcf} -Ou | \
		bcftools sort -Oz -o mutlisampleID_${params.runID}_mutect2_ROA_FAA_PASS.vcf.gz

		tabix -fp vcf mutlisampleID_${params.runID}_mutect2_ROA_FAA_PASS.vcf.gz

		"""	
}
