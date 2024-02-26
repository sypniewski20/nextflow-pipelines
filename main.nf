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
include { DELLY_CNV_CALL; DELLY_FILTER_VCF; DELLY_MERGE_SITES; DELLY_MERGED_SITES_CALL; DELLY_MERGED_FILTER } from "./modules/delly.nf"
include { MANTA_CNV_CALL; MANTA_FILTER_VCF } from "./modules/manta.nf"
include { ERDS_CNV_CALL; ERDS_FILTER_VCF; ERDS_FILTER_AND_MERGE_VCF } from "./modules/erds.nf"
include { CNVPYTOR_CALL; CNVPYTOR_FILTER_VCF; CNVPYTOR_MERGE_CALLS; CNVPYTOR_FILTER_MERGED } from "./modules/CNVpytor.nf"
include { WRITE_BAM_CHECKPOINT } from './modules/checkpoint.nf'
include { SV2_REGENOTYPE; SV2_FILTER_VCF; JOINT_SV2_REGENOTYPE; JOINT_SV2_FILTER_VCF } from './modules/sv2.nf'


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
			SAMBAMBA_MARK_DUPLICATES.out.ch_bam | set { ch_bam_filtered }
			WRITE_BAM_CHECKPOINT(SAMBAMBA_MARK_DUPLICATES.out.sample_checkpoint.collect())
		} else {
			FASTQ_TO_SAM(ch_fastp_results, fasta)
			BWA_SPARK_MAP_READS(FASTQ_TO_SAM.out, fasta, interval_list)		
			BQSR_SPARK(BWA_SPARK_MAP_READS.out, fasta, interval_list, snv_resource, cnv_resource, sv_resource)
			MARK_DUPLICATES_SPARK(BQSR_SPARK.out, interval_list)
			MARK_DUPLICATES_SPARK.out.ch_bam | set { ch_bam_filtered }
			WRITE_BAM_CHECKPOINT(MARK_DUPLICATES_SPARK.out.sample_checkpoint.collect())
		}
		WRITE_BAM_CHECKPOINT.out | set { ch_checkpoint }
		MOSDEPTH_WGS(ch_bam_filtered)
	emit:
		ch_bam_filtered
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

workflow paired_end_cnv_call {
	take:
		ch_samples_checkpoint
		fasta
		centromeres
	main:
		MANTA_CNV_CALL(ch_samples_checkpoint, fasta)
		MANTA_FILTER_VCF(MANTA_CNV_CALL.out)
		MANTA_FILTER_VCF.out | set { manta_output }
		if (params.cohort_mode == false) {
			DELLY_CNV_CALL(ch_samples_checkpoint, fasta, centromeres)
			DELLY_FILTER_VCF(DELLY_CNV_CALL.out)
			DELLY_FILTER_VCF.out | set { delly_output }
		} else {
			DELLY_CNV_CALL(ch_samples_checkpoint, fasta, centromeres)
			DELLY_MERGE_SITES(DELLY_CNV_CALL.out.collect())
			DELLY_MERGED_SITES_CALL(ch_samples_checkpoint, fasta, centromeres, DELLY_MERGE_SITES.out)
			DELLY_MERGED_FILTER(DELLY_MERGED_SITES_CALL.out.collect()) | set { delly_output }
		}
	emit:
		delly_output
		manta_output
}

workflow coverage_based_cnv_call {
	take:
		ch_samples_checkpoint
		fasta
		ch_snv_call
	main:
		if (params.cohort_mode == false) {
			CNVPYTOR_CALL(ch_samples_checkpoint)
			CNVPYTOR_FILTER_VCF(CNVPYTOR_CALL.out.vcf)
			CNVPYTOR_FILTER_VCF.out | set { cnvpytor_output }
			ch_samples_checkpoint
				.join(ch_snv_call)
				.set{erds_input}
			ERDS_CNV_CALL(erds_input, fasta)
			ERDS_FILTER_VCF(ERDS_CNV_CALL.out)
			ERDS_FILTER_VCF.out | set { erds_output }
		} else {
			CNVPYTOR_CALL(ch_samples_checkpoint)
			CNVPYTOR_MERGE_CALLS(CNVPYTOR_CALL.out.pytor.collect())
			CNVPYTOR_FILTER_MERGED(CNVPYTOR_MERGE_CALLS.out)
			CNVPYTOR_FILTER_MERGED.out | set { cnvpytor_output }
			ch_samples_checkpoint
				.join(ch_snv_call)
				.set{erds_input}
			ERDS_CNV_CALL(erds_input, fasta)
			ERDS_FILTER_AND_MERGE_VCF(ERDS_CNV_CALL.out.collect())
			ERDS_FILTER_AND_MERGE_VCF.out | set { erds_output }
		} 	
	emit:
		cnvpytor_output
		erds_output
}

workflow joint_cnv_regenotyping {
	take:
		ch_samples_checkpoint
		delly_output
		manta_output
		cnvpytor_output
		erds_output
		ch_snv_call
		ped_file
		fasta
	main:
		if (params.cohort_mode == false) {
			SV2_REGENOTYPE(ch_samples_checkpoint,
						delly_output, 
						manta_output,
						cnvpytor_output,
						erds_output,
						ch_snv_call, ped_file, fasta)
			SV2_FILTER_VCF(SV2_REGENOTYPE.out)
		} else {
			JOINT_SV2_REGENOTYPE(
						ch_samples_checkpoint.collect(),
						delly_output, 
						manta_output,
						cnvpytor_output,
						erds_output,
						ch_snv_call, ped_file, fasta)
			JOINT_SV2_FILTER_VCF(JOINT_SV2_REGENOTYPE.out)

		}

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
		bam_input = mapping_workflow.out.ch_bam_filtered
	} else {
		bam_checkpoint(params.bam_samplesheet) | set { bam_input }
	}
		snv_call(bam_input, 
				 params.fasta)
		paired_end_cnv_call(bam_input, 
							params.fasta,
							params.centromeres)
		coverage_based_cnv_call(bam_input, 
								params.fasta, 
								snv_call.out)
		joint_cnv_regenotyping(bam_input,
							   paired_end_cnv_call.out.delly_output,
							   paired_end_cnv_call.out.manta_output,
							   coverage_based_cnv_call.out.cnvpytor_output,
							   coverage_based_cnv_call.out.erds_output,
							   snv_call.out,
							   params.ped_file,
							   params.fasta)
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
