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
Intervals			:${params.interval_list}
SNV					:${params.snv_resource}
CNV					:${params.cnv_resource}
SV					:${params.sv_resource}
=================================================================================================================
"""
.stripIndent()

include { Read_samplesheet; Read_bam_checkpoint } from './modules/functions.nf'
include { FASTQC_PROCESSING; FASTP_PROCESSING; MOSDEPTH_WGS} from './modules/seqQC.nf'
include { BWA_MAP_READS; BASE_RECALIBRATOR; APPLY_BQSR; SAMBAMBA_MARK_DUPLICATES } from './modules/mapping.nf'
include { FASTQ_TO_SAM; BWA_SPARK_MAP_READS; BQSR_SPARK; MARK_DUPLICATES_SPARK } from './modules/spark_workflows.nf'
include { DEEP_VARIANT; FILTER_SNVS; FILTER_AND_MERGE_SNVS } from "./modules/deepvariant.nf"
include { MANTA_GERMLINE; MANTA_FILTER_VCF; MANTA_MERGE_VCF } from "./modules/manta.nf"
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
		snv_resource
		cnv_resource
		sv_resource
	main:
		if( params.run_spark == false ) {
			BWA_MAP_READS(ch_fastp_results, fasta)		
			BASE_RECALIBRATOR(BWA_MAP_READS.out, fasta, interval_list, snv_resource, cnv_resource, sv_resource)
			APPLY_BQSR(BASE_RECALIBRATOR.out, interval_list, fasta)
			SAMBAMBA_MARK_DUPLICATES(APPLY_BQSR.out)
			SAMBAMBA_MARK_DUPLICATES.out.ch_bam | set { ch_bam }
			WRITE_BAM_CHECKPOINT(SAMBAMBA_MARK_DUPLICATES.out.sample_checkpoint.collect())
		} else {
			FASTQ_TO_SAM(ch_fastp_results, fasta)
			BWA_SPARK_MAP_READS(FASTQ_TO_SAM.out, fasta, interval_list)		
			BQSR_SPARK(BWA_SPARK_MAP_READS.out, fasta, interval_list, snv_resource, cnv_resource, sv_resource)
			BQSR_SPARK.out.ch_bam | set { ch_bam }
			WRITE_BAM_CHECKPOINT(BQSR_SPARK.out.sample_checkpoint.collect())
		}
		WRITE_BAM_CHECKPOINT.out | set { ch_checkpoint }
		if( params.exome == false ) {
			MOSDEPTH_WGS(ch_bam_filtered)
		} else if( params.exome == false ) {
			MOSDEPTH_EXOME(ch_bam_filtered)
		}
	emit:
		ch_bam
		ch_checkpoint
}

workflow bam_checkpoint {
	take:
		ch_checkpoint
	main:
		Read_bam_checkpoint(ch_checkpoint)
	emit:
		ch_samples_checkpoint
}

workflow snv_call {
	take:
		ch_samples_checkpoint
		fasta
	main:
		if( params.cohort_mode == false ) {
			DEEP_VARIANT(ch_samples_checkpoint, fasta)
			FILTER_SNVS(DEEP_VARIANT.out)
			FILTER_SNVS.out | set { ch_snv_call }
		} else {
			DEEP_VARIANT(ch_samples_checkpoint, fasta)
			FILTER_AND_MERGE_SNVS(DEEP_VARIANT.out.collect(), fasta)
			FILTER_AND_MERGE_SNVS.out | set { ch_snv_call }
		}
	emit:
		ch_snv_call
}

workflow sv_call {
	take:
		ch_samples_checkpoint
		fasta
	main:
		if( params.exome == false ) {
			MANTA_GERMLINE(ch_samples_checkpoint, fasta)
			MANTA_FILTER_VCF(MANTA_GERMLINE.out)
		} else {
			MANTA_EXOME_GERMLINE(ch_samples_checkpoint, fasta)
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
		VEP(ch_snv_call, fasta)
		ANNOT_SV(ch_sv_output)
}

//Main workflow

workflow {
	main:
	if (params.bam_samplesheet.isEmpty()) {
		preprocessing_workflow(params.samplesheet)
		fastq_QC_workflow(preprocessing_workflow.out.ch_samples_initial)
		mapping_workflow(fastq_QC_workflow.out.ch_fastp_results, 
						 params.fasta, 
						 params.interval_list, 
						 params.snv_resource, 
						 params.cnv_resource, 
						 params.sv_resource)
		bam_input = mapping_workflow.out.ch_bam
	} else {
		bam_checkpoint(params.bam_samplesheet) | set { bam_input }
	}
		snv_call(bam_input, 
				 params.fasta)
		sv_call(bam_input, 
					params.fasta)
		annotation_workflow(snv_call.out, sv_call.out, params.fasta)
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
