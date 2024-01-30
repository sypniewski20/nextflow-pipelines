nextflow.enable.dsl = 2

println """\
=================================================================================================================
										Run details
=================================================================================================================
Samplesheet 		:${params.samplesheet}
Covarage filter		:Regions with > ${params.coverage_threshold}x coverage
Output				:${params.outfolder}
=================================================================================================================	
										References
=================================================================================================================
Fasta 				:${params.fasta}
Germline resource   :${params.germline_resource}
=================================================================================================================
"""
.stripIndent()

include { Read_samplesheet; Read_bam_checkpoint } from './modules/functions.nf'
include { CHECK_INTEGRITY; FASTQC_PROCESSING; FASTP_PROCESSING; BASE_RECALIBRATOR; APPLY_BQSR; MARK_DUPLICATES; MOSDEPTH_EXOME; COVERAGE_FILTER} from './modules/seqQC.nf'
include { CHECK_INTEGRITY; FASTQC_PROCESSING; FASTP_PROCESSING; MOSDEPTH_EXOME; COVERAGE_FILTER} from './modules/seqQC.nf'
include { BUILD_BAM_INDEX } from './modules/Mapping.nf'
include { FASTQ_TO_SAM; BWA_SPARK_MAP_READS; BASE_RECALIBRATOR_SPARK; APPLY_BQSR_SPARK; MARK_DUPLICATES_SPARK; WRITE_BAM_CHECKPOINT} from './modules/spark_workflows.nf'
include { MUTECT2_NORMAL; MUTECT2_PON; MUTECT2_ROA_MODEL; MUTECT2_PILEUP; MUTECT2_CONTAMINATION; Fun_to_contamination_prep; MUTECT2_ROA_FILTER; MUTECT2_AS_FILTER_HEADER; FILTER_AND_MERGE_VCF } from './modules/Mutect2_calling.nf'
include { SIG_INPUT_PREP; SIG_EXTRACTOR } from './modules/signatures.nf'


workflow preprocessing_workflow {
	take:
		samplesheet
	main:
		Read_samplesheet(samplesheet)
	emit:
		ch_samples_initial
}

workflow fastq_QC_workflow {
	take:
		ch_samples_initial
	main:
		CHECK_INTEGRITY(ch_samples_initial)
		FASTQC_PROCESSING(CHECK_INTEGRITY.out) | set { fastqc_log }
		FASTP_PROCESSING(CHECK_INTEGRITY.out)
		FASTP_PROCESSING.out.fastq_filtered | set { ch_fastp_results }
	emit:
		fastqc_log
		ch_fastp_results
}

workflow mapping_workflow {
	take:
		ch_fastp_results
		fasta
		interval_list
		interval_padding
		germline_resource
		exome_bed
	main:
		FASTQ_TO_SAM(ch_fastp_results, fasta)
		BWA_SPARK_MAP_READS(FASTQ_TO_SAM.out, interval_list, fasta, interval_padding)		
		BASE_RECALIBRATOR_SPARK(BWA_SPARK_MAP_READS.out, interval_list, fasta, germline_resource)
		BWA_SPARK_MAP_READS.out
			.combine(BASE_RECALIBRATOR_SPARK.out) | set { bqsr_input }
		APPLY_BQSR_SPARK(bqsr_input, interval_list, fasta)
		MARK_DUPLICATES_SPARK(APPLY_BQSR_SPARK.out, interval_list)
		MARK_DUPLICATES_SPARK.out.ch_bam | set { ch_bam_filtered }
		WRITE_BAM_CHECKPOINT(MARK_DUPLICATES_SPARK.out.sample_checkpoint.collect())
		MOSDEPTH_EXOME(ch_bam_filtered, exome_bed)
		COVERAGE_FILTER(MOSDEPTH_EXOME.out.collect())
		COVERAGE_FILTER.out.coverage_bed | set { ch_coverage_bed }
	emit:
		ch_bam_filtered
		ch_coverage_bed
}

workflow bam_checkpoint {
	take:
		bam_checkpoint_sheet
	main:
		Read_samplesheet(bam_checkpoint_sheet)
	emit:
		ch_samples_initial
}


workflow pon_workflow {
	take:
		ch_bam_filtered
		fasta
		ch_coverage_bed
		interval_padding
	main:
		MUTECT2_NORMAL(ch_bam_filtered, fasta, ch_coverage_bed, interval_padding)
		MUTECT2_PON(MUTECT2_NORMAL.out.collect(), fasta, ch_coverage_bed, interval_padding) | set { ch_pon }
	emit:
		ch_pon
}

workflow snv_somatic_call_workflow {
	take:
		ch_bam_filtered
		fasta
		germline_resource
		ch_coverage_bed
		ch_pon
		interval_padding
	main:
		MUTECT2_ROA_MODEL(ch_bam_filtered, fasta, germline_resource, ch_coverage_bed, ch_pon, interval_padding)
		MUTECT2_PILEUP(ch_bam_filtered, germline_resource, ch_coverage_bed, interval_padding)
		MUTECT2_CONTAMINATION(MUTECT2_PILEUP.out)
		Fun_to_contamination_prep(MUTECT2_ROA_MODEL.out.roa_model, MUTECT2_CONTAMINATION.out)
		MUTECT2_ROA_FILTER(ch_cont_input, fasta)
		MUTECT2_AS_FILTER_HEADER(MUTECT2_ROA_FILTER.out, fasta)
		MUTECT2_AS_FILTER_HEADER.out.vcf.collect() | set { ch_sigs_input }
		FILTER_AND_MERGE_VCF(MUTECT2_AS_FILTER_HEADER.out.vcf.collect(), MUTECT2_AS_FILTER_HEADER.out.tbi.collect()) | set { ch_somatic_output }
	emit:
		ch_somatic_output
		ch_sigs_input
}

workflow signatures_workflow {
	take:
		ch_sigs_input
	main:
		SIG_INPUT_PREP(ch_sigs_input)
		SIG_EXTRACTOR(SIG_INPUT_PREP.out) | set { sig_output }
	emit:
		sig_output
}
	
//Main workflow
workflow {
	main:
	if (params.bam_checkpoint.isEmpty()) {
		preprocessing_workflow(params.samplesheet)
		fastq_QC_workflow(preprocessing_workflow.out.ch_samples_initial)
		mapping_workflow(fastq_QC_workflow.out.ch_fastp_results, params.fasta, params.interval_list, params.interval_padding, params.germline_resource, params.exome_bed)
		pon_workflow(mapping_workflow.out.ch_bam_filtered, params.fasta, mapping_workflow.out.ch_coverage_bed, params.interval_padding)
		snv_somatic_call_workflow(mapping_workflow.out.ch_bam_filtered, params.fasta, params.germline_resource, mapping_workflow.out.ch_coverage_bed, pon_workflow.out.ch_pon, params.interval_padding)
	} else {
		bam_checkpoint(params.bam_checkpoint)
		mapping_workflow(bam_checkpoint.out.ch_samples_initial, params.fasta, params.interval_list, params.interval_padding, params.germline_resource, params.exome_bed)
		pon_workflow(bam_checkpoint.out.ch_samples_initial, params.fasta, mapping_workflow.out.ch_coverage_bed)
		snv_somatic_call_workflow(bam_checkpoint.out.ch_samples_initial, params.fasta, params.germline_resource, mapping_workflow.out.ch_coverage_bed, pon_workflow.out.ch_pon)
		signatures_workflow(snv_somatic_call_workflow.out.ch_sigs_input)
	}

		signatures_workflow(snv_somatic_call_workflow.out.ch_sigs_input)
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
