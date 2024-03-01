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
include { CHECK_INTEGRITY; FASTQC_PROCESSING; FASTP_PROCESSING; MOSDEPTH_WGS} from './modules/seqQC.nf'
include { BWA_MAP_READS; BASE_RECALIBRATOR; APPLY_BQSR; SAMBAMBA_MARK_DUPLICATES } from './modules/mapping.nf'
include { FASTQ_TO_SAM; BWA_SPARK_MAP_READS; BQSR_SPARK; MARK_DUPLICATES_SPARK } from './modules/spark_workflows.nf'
include { DEEP_VARIANT; FILTER_SNVS; FILTER_AND_MERGE_SNVS } from "./modules/deepvariant.nf"
include { MANTA_CNV_CALL; MANTA_FILTER_VCF } from "./modules/manta.nf"
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
		MOSDEPTH_WGS(ch_bam)
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
			FILTER_SNVS(DEEP_VARIANT.out, fasta)
			FILTER_SNVS.out | set { ch_snv_call }
		} else {
			DEEP_VARIANT(ch_samples_checkpoint, fasta)
			FILTER_AND_MERGE_SNVS(DEEP_VARIANT.out.collect(), fasta)
			FILTER_AND_MERGE_SNVS.out | set { ch_snv_call }
		}
	emit:
		ch_snv_call
}

workflow manta_call {
	take:
		ch_samples_checkpoint
		fasta
	main:
		MANTA_CNV_CALL(ch_samples_checkpoint, fasta)
		MANTA_FILTER_VCF(MANTA_CNV_CALL.out)
		MANTA_FILTER_VCF.out | set { manta_output }
	emit:
		manta_output
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
		manta_call(bam_input, 
					params.fasta)
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
