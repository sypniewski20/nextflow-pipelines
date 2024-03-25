nextflow.enable.dsl = 2

println """\
=================================================================================================================
										Run details: 
=================================================================================================================
Run ID              :${params.runID}
Samplesheet 		:${params.samplesheet}
Output				:${params.outfolder}
=================================================================================================================	
										References
=================================================================================================================
Fasta 				:${params.fasta}
=================================================================================================================
"""
.stripIndent()

include { Read_samplesheet; Read_bam_checkpoint } from './modules/functions.nf'
include { FASTQC_PROCESSING; FASTP_PROCESSING; MOSDEPTH_WGS; MOSDEPTH_EXOME} from './modules/seqQC.nf'
include { BWA_MAP_READS } from './modules/mapping.nf'
include { BASE_RECALIBRATOR; APPLY_BQSR; BQSR_SPARK } from './modules/bqsr.nf'
include { DEEP_VARIANT; FILTER_SNVS; FILTER_AND_MERGE_SNVS } from "./modules/deepvariant.nf"
include { MANTA_GERMLINE; MANTA_EXOME_GERMLINE; MANTA_FILTER_VCF; MANTA_MERGE_VCF } from "./modules/manta.nf"
include { VEP; ANNOT_SV } from "./modules/annotations.nf"
include { WRITE_BAM_CHECKPOINT } from './modules/checkpoint.nf'


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
		FASTQC_PROCESSING(ch_samples_initial) | set { fastqc_log }
		FASTP_PROCESSING(ch_samples_initial)
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
	main:
		if (params.bwa == 1) {
			BWA_MAP_READS(ch_fastp_results, fasta) | set { ch_bam }	
		} else if (params.bwa == 2) {
			BWAMEM2_MAP_READS(ch_fastp_results, fasta) | set { ch_bam }
		}
		
		if( params.exome == false ) {
			MOSDEPTH_WGS(ch_bam)
		} else if( params.exome == true ) {
			MOSDEPTH_EXOME(ch_bam)
		}
	emit:
		ch_bam
}

workflow bqsr_mapping_workflow {
	take:
		ch_fastp_results
		fasta
		interval_list
		snv_resource
	main:
		if (params.bwa == 1) {
			BWA_MAP_READS(ch_fastp_results, fasta) | set { ch_bam }	
		} else if (params.bwa == 2) {
			BWAMEM2_MAP_READS(ch_fastp_results, fasta) | set { ch_bam }
		}

		BASE_RECALIBRATOR(ch_bam.out, fasta, interval_list, snv_resource)
		APPLY_BQSR(BASE_RECALIBRATOR.out, interval_list, fasta)

		if( params.exome == false ) {
			MOSDEPTH_WGS(ch_bam)
		} else if( params.exome == true ) {
			MOSDEPTH_EXOME(ch_bam)
		}
	emit:
		ch_bam
}

workflow read_bam {
	take:
		ch_checkpoint
	main:
		Read_bam_checkpoint(ch_checkpoint)
	emit:
		ch_bam
}

workflow snv_call {
	take:
		ch_bam
		fasta
	main:
		if( params.cohort_mode == false ) {
			DEEP_VARIANT(ch_bam, fasta)
			FILTER_SNVS(DEEP_VARIANT.out)
			FILTER_SNVS.out | set { ch_snv_call }
		} else {
			DEEP_VARIANT(ch_bam, fasta)
			FILTER_AND_MERGE_SNVS(DEEP_VARIANT.out.collect(), fasta)
			FILTER_AND_MERGE_SNVS.out | set { ch_snv_call }
		}
	emit:
		ch_snv_call
}

workflow sv_call {
	take:
		ch_bam
		fasta
	main:
		if( params.exome == false ) {
			MANTA_GERMLINE(ch_bam, fasta)
			MANTA_FILTER_VCF(MANTA_GERMLINE.out)
		} else {
			MANTA_EXOME_GERMLINE(ch_bam, fasta)
			MANTA_FILTER_VCF(MANTA_EXOME_GERMLINE.out)
		}
		if( params.cohort_mode == false) {
			MANTA_FILTER_VCF.out | set { ch_sv_output }
		} else {
			MANTA_FILTER_VCF.out
				.map{ sample, vcf, tbi -> vcf }
				.collect()
				.set{ vcf_input }
			MANTA_FILTER_VCF.out
				.map{ sample, vcf, tbi -> tbi }
				.collect()
				.set{ tbi_input }
			MANTA_MERGE_VCF(vcf_input, tbi_input) | set { ch_sv_output }
		}
	emit:
		ch_sv_output
}

workflow annotation_workflow {
	take:
		ch_snv_call
		ch_sv_output
		fasta
	main:
		// VEP(ch_snv_call, fasta)
		ANNOT_SV(ch_sv_output)
}

//Main workflow

workflow {
	main:
	if (params.mode == "mapping" ) {

		preprocessing_workflow(params.samplesheet)

		fastq_QC_workflow(preprocessing_workflow.out.ch_samples_initial)

		if (params.bqsr == false) {
			mapping_workflow(fastq_QC_workflow.out.ch_fastp_results, 
						params.fasta, 
						params.interval_list)

		} else {
			bqsr_mapping_workflow(fastq_QC_workflow.out.ch_fastp_results, 
						params.fasta, 
						params.interval_list, 
						params.snv_resource)

		}
		
	} else if(params.mode == "calling" ) {

		read_bam(params.bam_samplesheet) | set { bam_input }

		snv_call(bam_input, 
				 params.fasta)

		sv_call(bam_input, 
				params.fasta)

		annotation_workflow(snv_call.out, sv_call.out, params.fasta)

	} else if(params.mode == "both" ) {

		preprocessing_workflow(params.samplesheet)
		fastq_QC_workflow(preprocessing_workflow.out.ch_samples_initial)

		if (params.bqsr == false) {
			mapping_workflow(fastq_QC_workflow.out.ch_fastp_results, 
						params.fasta, 
						params.interval_list) | set { ch_bam }

		} else {
			bqsr_mapping_workflow(fastq_QC_workflow.out.ch_fastp_results, 
						params.fasta, 
						params.interval_list, 
						params.snv_resource) | set { ch_bam }

		}

		snv_call(ch_bam, 
				 params.fasta)

		sv_call(ch_bam, 
					params.fasta)

		annotation_workflow(snv_call.out, sv_call.out, params.fasta)
	}
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
