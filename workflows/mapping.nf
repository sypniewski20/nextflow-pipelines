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

include { Read_samplesheet } from '../modules/functions.nf'
include { FASTQC_PROCESSING; FASTP_PROCESSING; MOSDEPTH_WGS; MOSDEPTH_EXOME} from '../modules/seqQC.nf'
include { BWA_MAP_READS; BWAMEM2_MAP_READS } from '../modules/mapping.nf'
include { WRITE_BAM_CHECKPOINT } from '../modules/checkpoint.nf'


workflow preprocessing_workflow {
	take:
		samplesheet
	main:
		Read_samplesheet(samplesheet)
	emit:
		ch_samples_initial
}

workflow read_bam {
	take:
		ch_checkpoint
	main:
		Read_bam_checkpoint(ch_checkpoint)
	emit:
		ch_bam
}

workflow bam2fastq {
	take:
		ch_bam
	main:
		UBAM(ch_bam)
		BAM2FASTQ(UBAM.out) 
		BAM2FASTQ.out | set { ch_samples_initial }
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


//Main workflow

workflow {
	main:
		if(params.remap == true) {
			read_bam(params.samplesheet)
			bam2fastq(read_bam.out)
			fastq_QC_workflow(bam2fastq.out)
		} else {
			preprocessing_workflow(params.samplesheet)
			fastq_QC_workflow(preprocessing_workflow.out.ch_samples_initial)
		}

        mapping_workflow(fastq_QC_workflow.out.ch_fastp_results, 
                    params.fasta)
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
