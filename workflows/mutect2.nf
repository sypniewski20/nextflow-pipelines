nextflow.enable.dsl = 2

println """\
=================================================================================================================
										Run details
=================================================================================================================
Samplesheet 		:${params.samplesheet}
Output				:${params.outfolder}
=================================================================================================================	
										References
=================================================================================================================
Fasta 				:${params.fasta}
SNV resource   :${params.snv_resource}
=================================================================================================================
"""
.stripIndent()

include { Read_samplesheet; Read_bam_checkpoint } from './modules/functions.nf'
include { CHECK_INTEGRITY; FASTQC_PROCESSING; FASTP_PROCESSING; MOSDEPTH_EXOME; COVERAGE_FILTER} from './modules/seqQC.nf'
include { BWA_MAP_READS; BASE_RECALIBRATOR; APPLY_BQSR; SAMBAMBA_MARK_DUPLICATES } from './modules/mapping.nf'
include { WRITE_BAM_CHECKPOINT } from './modules/checkpoint.nf'
include { MUTECT2_NORMAL; MUTECT2_PON; MUTECT2_ROA_MODEL; MUTECT2_PILEUP; MUTECT2_CONTAMINATION; Fun_to_contamination_prep; MUTECT2_ROA_FILTER; MUTECT2_AS_FILTER_HEADER; MERGE_VCF } from './modules/Mutect2_calling.nf'
include { MANTA_EXOME_TUMOR_ONLY; MANTA_FILTER_VCF; MANTA_MERGE_VCF } from "./modules/manta.nf"
include { SIG_INPUT_PREP; SIG_DIR_PREP; SIG_EXTRACTOR } from './modules/signatures.nf'
include { ANNOT_SV } from './modules/annotsv.nf'


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
		FASTP_PROCESSING.out.fastp_log | set { fastp_log }
	emit:
		fastqc_log
		fastp_log
		ch_fastp_results
}

workflow mapping_workflow {
	take:
		ch_fastp_results
		fasta
		interval_list
		interval_padding
		snv_resource
	main:
		BWA_MAP_READS(ch_fastp_results, fasta)		
		BASE_RECALIBRATOR(BWA_MAP_READS.out, fasta, interval_list, interval_padding, snv_resource)
		APPLY_BQSR(BASE_RECALIBRATOR.out, fasta, interval_padding)
		SAMBAMBA_MARK_DUPLICATES(APPLY_BQSR.out)
		SAMBAMBA_MARK_DUPLICATES.out.ch_bam | set { ch_bam_filtered }
		WRITE_BAM_CHECKPOINT(SAMBAMBA_MARK_DUPLICATES.out.sample_checkpoint.collect())
		WRITE_BAM_CHECKPOINT.out | set { ch_checkpoint }
	emit:
		ch_bam_filtered
		ch_checkpoint
}

workflow mosdepth_workflow {
	take:
		ch_bam_filtered
		interval_list
	main:
		MOSDEPTH_EXOME(ch_bam_filtered, interval_list)
	// 	COVERAGE_FILTER(MOSDEPTH_EXOME.out.collect())
	// 	COVERAGE_FILTER.out.coverage_bed | set { ch_coverage_bed }
	// emit:
	// 	ch_coverage_bed
}

workflow bam_checkpoint {
	take:
		ch_checkpoint
	main:
		Read_bam_checkpoint(ch_checkpoint)
	emit:
		ch_samples_checkpoint
}

workflow manta_somatic_only {
	take:
		ch_samples_checkpoint
		fasta
	main:
		MANTA_EXOME_TUMOR_ONLY(ch_samples_checkpoint, fasta)
		MANTA_FILTER_VCF(MANTA_EXOME_TUMOR_ONLY.out)
		MANTA_FILTER_VCF.out | set { manta_output }
		MANTA_FILTER_VCF.out
			.map{ sample, vcf, tbi -> [vcf] }
			.collect()
			.set{ vcf_input }
		MANTA_FILTER_VCF.out
			.map{ sample, vcf, tbi -> [tbi] }
			.collect()
			.set{ tbi_input }
		MANTA_MERGE_VCF(vcf_input, tbi_input)
		ANNOT_SV(MANTA_MERGE_VCF.out)
	emit:
		manta_output
}

workflow pon_workflow {
	take:
		ch_samples_checkpoint
		fasta
		ch_coverage_bed
		interval_padding
	main:
		MUTECT2_NORMAL(ch_samples_checkpoint, fasta, ch_coverage_bed, interval_padding)
		MUTECT2_PON(MUTECT2_NORMAL.out.collect(), fasta, ch_coverage_bed, interval_padding) | set { ch_pon }
	emit:
		ch_pon
}

workflow snv_somatic_call_workflow {
	take:
		ch_samples_checkpoint
		fasta
		snv_resource
		ch_coverage_bed
		ch_pon
		interval_padding
	main:
		MUTECT2_ROA_MODEL(ch_samples_checkpoint, fasta, snv_resource, ch_coverage_bed, ch_pon, interval_padding)
		MUTECT2_PILEUP(ch_samples_checkpoint, snv_resource, ch_coverage_bed, interval_padding)
		MUTECT2_CONTAMINATION(MUTECT2_PILEUP.out)
		Fun_to_contamination_prep(MUTECT2_ROA_MODEL.out, MUTECT2_CONTAMINATION.out) | set{ ch_cont_input }
		MUTECT2_ROA_FILTER(ch_cont_input, fasta)
		MUTECT2_AS_FILTER_HEADER(MUTECT2_ROA_FILTER.out, fasta)
		MERGE_VCF(MUTECT2_AS_FILTER_HEADER.out.vcf.collect(), MUTECT2_AS_FILTER_HEADER.out.tbi.collect()) | set { ch_somatic_output }
	emit:
		ch_somatic_output
}

workflow signatures_workflow {
	take:
		vcf
		fasta
	main:
		SIG_INPUT_PREP(vcf)
		SIG_DIR_PREP(SIG_INPUT_PREP.out.collect())
		SIG_EXTRACTOR(SIG_INPUT_PREP.out.collect()) | set { sig_output }
	emit:
		sig_output
}
	
//Main workflow
workflow {
	main:
	if (params.bam_samplesheet.isEmpty() ) {
		preprocessing_workflow(params.samplesheet)
		fastq_QC_workflow(preprocessing_workflow.out.ch_samples_initial)
		mapping_workflow(fastq_QC_workflow.out.ch_fastp_results, 
						 params.fasta, 
						 params.interval_list, 
						 params.interval_padding, 
						 params.snv_resource)
		mapping_workflow.out.ch_bam_filtered | set { bam_input }
		mosdepth_workflow(bam_input, params.interval_list)
	} else {
		bam_checkpoint(params.bam_samplesheet) | set { bam_input }
	}
		manta_somatic_only(bam_input, params.fasta)
		pon_workflow(bam_input, 
					 params.fasta,
					 params.interval_list,  
					 params.interval_padding)
		snv_somatic_call_workflow(bam_input, 
								  params.fasta, 
								  params.snv_resource, 
								  params.interval_list, 
								  pon_workflow.out.ch_pon, 
								  params.interval_padding)
		signatures_workflow(snv_somatic_call_workflow.out, params.fasta)
}


// workflow {
// 	main:
// 	if (params.bam_samplesheet.isEmpty() ) {
// 		preprocessing_workflow(params.samplesheet)
// 		fastq_QC_workflow(preprocessing_workflow.out.ch_samples_initial)
// 		mapping_workflow(fastq_QC_workflow.out.ch_fastp_results, 
// 						 params.fasta, 
// 						 params.interval_list, 
// 						 params.interval_padding, 
// 						 params.snv_resource)
// 		mapping_workflow.out.ch_bam_filtered | set { bam_input }
// 		mosdepth_workflow(bam_input, params.interval_list)
// 	} else {
// 		bam_checkpoint(params.bam_samplesheet) | set { bam_input }
// 		mosdepth_workflow(bam_input, params.interval_list)
// 	}
// 		pon_workflow(bam_input, 
// 					 params.fasta, 
// 					 mosdepth_workflow.out.ch_coverage_bed, 
// 					 params.interval_padding)
// 		snv_somatic_call_workflow(bam_input, 
// 								  params.fasta, 
// 								  params.snv_resource, 
// 								  mosdepth_workflow.out.ch_coverage_bed, 
// 								  pon_workflow.out.ch_pon, 
// 								  params.interval_padding)
// 		signatures_workflow(snv_somatic_call_workflow.out.ch_sigs_input)
// }

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
