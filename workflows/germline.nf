nextflow.enable.dsl = 2

println """\
=================================================================================================================
										Run details: 
=================================================================================================================
Run ID              :${params.runID}
Output				:${params.outfolder}/${params.runID}
=================================================================================================================	
										References
=================================================================================================================
Fasta 				:${params.fasta}
=================================================================================================================
"""
.stripIndent()

include { Read_samplesheet; Read_bam_checkpoint } from '../modules/functions.nf'
include { FASTQC_PROCESSING; FASTP_PROCESSING; MOSDEPTH_WGS; MOSDEPTH_EXOME} from '../modules/seqQC.nf'
include { BWA_MAP_READS; BWAMEM2_MAP_READS } from '../modules/mapping.nf'
include { BASE_RECALIBRATOR; APPLY_BQSR; BQSR_SPARK } from '../modules/bqsr.nf'
include { DEEP_VARIANT_WGS; DEEP_VARIANT_WES; FILTER_SNVS } from "../modules/deepvariant.nf"
include { MANTA_GERMLINE; MANTA_EXOME_GERMLINE; MANTA_FILTER_VCF } from "../modules/manta.nf"
include { SMOOVE; SMOOVE_ANNOTATE } from "../modules/smoove.nf"
include { DELLY_SV_CALL; DELLY_CNV_CALL; FILTER_DELLY; DELLY_MERGE_SITES; DELLY_MERGED_SITES_CALL } from "../modules/delly.nf"
include { VEP_SNV; VEP_SV; ANNOT_SV; SURVIVOR } from "../modules/annotations.nf"
include { WRITE_BAM_CHECKPOINT } from '../modules/checkpoint.nf'


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
	main:
		fasta_fai = file(fasta+'.fai')
		fasta_sa = file(fasta+'.sa')
		fasta_bwt = file(fasta+'.bwt')
		fasta_ann = file(fasta+'.ann')
		fasta_amb = file(fasta+'.amb')
		fasta_pac = file(fasta+'.pac')
		fasta0123 = file(fasta+'.0123')
		fasta_2bit = file(fasta+'.bwt.2bit.64')

		if (params.bwa == 1) {
			BWA_MAP_READS(ch_fastp_results,
						  fasta,
						  fasta_fai,
						  fasta_sa,
						  fasta_bwt,
						  fasta_ann,
						  fasta_amb,
						  fasta_pac) | set { ch_bam }	
		} else if (params.bwa == 2) {

			BWAMEM2_MAP_READS(ch_fastp_results, 
							  fasta,
							  fasta_fai,
							  fasta0123,
							  fasta_2bit,
							  fasta_sa,
							  fasta_bwt,
							  fasta_ann,
							  fasta_amb,
							  fasta_pac) | set { ch_bam }
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
		fasta_fai = file(fasta+'.fai')
		fasta_sa = file(fasta+'.sa')
		fasta_bwt = file(fasta+'.bwt')
		fasta_ann = file(fasta+'.ann')
		fasta_amb = file(fasta+'.amb')
		fasta_pac = file(fasta+'.pac')
		fasta_dict = file(fasta-'.fa'+'.dict')
		snv_tbi = file(snv_resource+'.tbi')
		
		if (params.bwa == 1) {
			BWA_MAP_READS(ch_fastp_results,
						  fasta,
						  fasta_fai,
						  fasta_sa,
						  fasta_bwt,
						  fasta_ann,
						  fasta_amb,
						  fasta_pac) | set { ch_bam }	
		} else if (params.bwa == 2) {
			BWAMEM2_MAP_READS(ch_fastp_results, fasta) | set { ch_bam }
		}

		BASE_RECALIBRATOR(ch_bam,
						  fasta,
						  fasta_fai,
						  fasta_dict, 
						  interval_list, 
						  snv_resource,
						  snv_tbi)
						  
		APPLY_BQSR(BASE_RECALIBRATOR.out, 
						  interval_list, 
						  fasta,
						  fasta_fai,
						  fasta_dict)

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
		contigs_bed
	main:
		fasta_fai = file(fasta+'.fai')
		if( params.exome == false ) {
			DEEP_VARIANT_WGS(ch_bam, 
							fasta,
							fasta_fai,
							contigs_bed)
			FILTER_SNVS(DEEP_VARIANT_WGS.out, fasta)
		} else {
			DEEP_VARIANT_WES(ch_bam, 
							fasta,
							fasta_fai,
							contigs_bed)
			FILTER_SNVS(DEEP_VARIANT_WES.out, fasta)
		}
	VEP_SNV(FILTER_SNVS.out, fasta, fasta_fai)
}

workflow sv_manta {
	take:
		ch_bam
		fasta
		contigs_bed
	main:
		fasta_fai = file(fasta+'.fai')
		contigs_tbi = file(contigs_bed+".tbi")

		if( params.exome == false ) {
			MANTA_GERMLINE(ch_bam, 
						   fasta, 
						   fasta_fai, 
						   contigs_bed, 
						   contigs_tbi)
			MANTA_FILTER_VCF(MANTA_GERMLINE.out)
		} else {
			MANTA_EXOME_GERMLINE(ch_bam, 
								 fasta, 
								 fasta_fai,
								 contigs_bed,
								 contigs_tbi)
			MANTA_FILTER_VCF(MANTA_EXOME_GERMLINE.out)
		}	
		MANTA_FILTER_VCF.out | set { manta_sv }
	emit:
		manta_sv
}

workflow sv_delly {
	take:
		ch_bam
		fasta
		contigs_bed
		delly_map
	main:
		fasta_fai = file(fasta+'.fai')
		DELLY_SV_CALL(ch_bam, 
						fasta, 
						fasta_fai, 
						contigs_bed)
		FILTER_DELLY(DELLY_SV_CALL.out) | set {delly_sv}
	emit:
		delly_sv
}

workflow sv_smoove {
	take:
		ch_bam
		fasta
	main:
		fasta_fai = file(fasta+'.fai')
		SMOOVE(ch_bam, fasta, fasta_fai) | set { smoove_sv }
	emit:
		smoove_sv
}

workflow survivor {
	take:
		manta_sv
		delly_sv
		smoove_sv
		fasta
	main:
		fasta_fai = file(fasta+'.fai')
		SURVIVOR(manta_sv, delly_sv, smoove_sv)
		ANNOT_SV(SURVIVOR.out)

}

//Main workflow

workflow {
	main:
	if (params.mode == "mapping" ) {

		preprocessing_workflow(params.samplesheet)

		fastq_QC_workflow(preprocessing_workflow.out.ch_samples_initial)

		mapping_workflow(fastq_QC_workflow.out.ch_fastp_results, 
					params.fasta)

	} else if(params.mode == "calling" ) {

		read_bam(params.bam_samplesheet) | set { ch_bam }

		snv_call(read_bam.out.ch_bam, 
				 params.fasta)

		sv_manta(ch_bam, 
				params.fasta,
				params.contigs)

		sv_delly(ch_bam, 
				params.fasta,
				params.centromeres,
				params.delly_map)

		sv_smoove(ch_bam,
				  params.fasta)

		survivor(sv_manta.out,
				 sv_delly.out,
				 sv_smoove.out)


	} else if(params.mode == "all" ) {

		preprocessing_workflow(params.samplesheet)
		fastq_QC_workflow(preprocessing_workflow.out.ch_samples_initial)

		mapping_workflow(fastq_QC_workflow.out.ch_fastp_results, 
					params.fasta) | set { ch_bam }

		snv_call(ch_bam, 
				 params.fasta,
				 params.contigs_bed)

		sv_manta(ch_bam, 
				params.fasta,
				params.contigs_bed)

		sv_delly(ch_bam, 
				params.fasta,
				params.contigs_bed,
				params.delly_map)

		sv_smoove(ch_bam,
				  params.fasta)

		survivor(sv_manta.out,
				 sv_delly.out,
				 sv_smoove.out,
				 params.fasta)
	}
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
