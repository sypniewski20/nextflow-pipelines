nextflow.enable.dsl = 2

println """\
=================================================================================================================
                                    Run details: 
=================================================================================================================
Run ID:                             ${params.runID}
Output:                             ${params.outfolder}/${params.runID}
Joint call: 						${params.joint_call}
Exome:                              ${params.exome}
Manta:                              ${params.manta}
Delly:                              ${params.delly}
Smoove:                             ${params.smoove}
GRIDSS:                             ${params.gridss}
VirusBreakend:						${params.virusbreakend}
=================================================================================================================	
                                    References
=================================================================================================================
Fasta: 				                ${params.fasta}
Exclude BED:			            ${params.exclude_bed}
=================================================================================================================
"""
.stripIndent()

include { Read_bam_checkpoint } from '../modules/functions.nf'
include { MANTA_GERMLINE; MANTA_EXOME_GERMLINE; MANTA_WGS_JOINT; MANTA_WES_JOINT; MANTA_FILTER_VCF } from "../modules/manta.nf"
include { SMOOVE; SMOOVE_JOINT; SMOOVE_ANNOTATE; SMOOVE_FILTER_VCF } from "../modules/smoove.nf"
include { DELLY_SV_CALL; DELLY_CNV_CALL } from "../modules/delly.nf"
include { GRIDSS; GRIDSS_JOINT; VIRUS_BREAKEND; GRIDSS_FILTER_VCF } from "../modules/gridss.nf"
include { VEP_SV; ANNOT_SV; SURVIVOR } from "../modules/annotations.nf"

workflow read_bam {
	take:
		ch_checkpoint
	main:
		Read_bam_checkpoint(ch_checkpoint)
	emit:
		ch_bam
}


workflow manta {
	take:
		ch_bam
		fasta
	main:
		fasta_fai = file(fasta+'.fai')

		if( params.exome == false && params.joint_call == true) {

			ch_bam.map { sample, bam, bai -> bam }
				  .collect()
				  .set{ joint_bam }

			ch_bam.map { sample, bam, bai -> bai }
				  .collect()
				  .set{ joint_bai }

			MANTA_WGS_JOINT(joint_bam, 
						   joint_bai,
						   fasta, 
						   fasta_fai) | set { manta_out }

		} else if ( params.exome == true && params.joint_call == true) {

			ch_bam.map { sample, bam, bai -> bam }
				  .collect()
				  .set{ joint_bam }

			ch_bam.map { sample, bam, bai -> bai }
				  .collect()
				  .set{ joint_bai }

			MANTA_WES_JOINT(joint_bam, 
						   joint_bai,
						   fasta, 
						   fasta_fai) | set { manta_out }

		} else if ( params.exome == true && params.joint_call == false) {

			MANTA_EXOME_GERMLINE(ch_bam, 
						   fasta, 
						   fasta_fai) | set { manta_out }

		} else if ( params.exome == false && params.joint_call == false) {

			MANTA_GERMLINE(ch_bam, 
								 fasta, 
								 fasta_fai) | set { manta_out }
		}	

		MANTA_FILTER_VCF(manta_out) | set { manta_sv }
		VEP_SV(manta_sv, fasta, fasta_fai)

	emit:
		manta_sv
}

workflow delly {
	take:
		ch_bam
		fasta
		delly_map
	main:
		fasta_fai = file(fasta+'.fai')
		DELLY_SV_CALL(ch_bam, 
						fasta, 
						fasta_fai)
		DELLY_SV_CALL.out | set { delly_sv }

	emit:
		delly_sv
}

workflow smoove {
	take:
		ch_bam
		fasta
	main:
		fasta_fai = file(fasta+'.fai')

		if (params.joint_call == true) {

			ch_bam.map { sample, bam, bai -> bam }
				  .collect()
				  .set{ joint_bam }

			ch_bam.map { sample, bam, bai -> bai }
				  .collect()
				  .set{ joint_bai }

			SMOOVE_JOINT(joint_bam, 
				   joint_bai,
				   fasta, 
				   fasta_fai) | set { smoove_output }

		} else {

			SMOOVE(ch_bam, fasta, fasta_fai)  | set { smoove_output }

		}

		SMOOVE_FILTER_VCF(smoove_output) | set { smoove_sv }

		VEP_SV(smoove_sv, fasta, fasta_fai)

	emit:
		smoove_sv
}

workflow gridss {
	take:
		ch_bam
		fasta
	main:
		fasta_sa = file(fasta+'.sa')
		fasta_bwt = file(fasta+'.bwt')
		fasta_ann = file(fasta+'.ann')
		fasta_amb = file(fasta+'.amb')
		fasta_pac = file(fasta+'.pac')
		fasta_dict = file(fasta+'.dict')
		fasta_fai = file(fasta+'.fai')

		if (params.joint_call == true) {
			
			ch_bam.map { sample, bam, bai -> bam }
				  .collect()
				  .set{ joint_bam }

			ch_bam.map { sample, bam, bai -> bai }
				  .collect()
				  .set{ joint_bai }

			GRIDSS_JOINT_PREPROCESS(joint_bam,
				joint_bai, 						  
				fasta,
				fasta_fai,
				fasta_dict,
				fasta_sa,
				fasta_bwt,
				fasta_ann,
				fasta_amb,
				fasta_pac) | set { gridss_output }

		} else {

			GRIDSS(ch_bam, 						  
               fasta,
			   fasta_fai,
               fasta_dict,
			   fasta_sa,
			   fasta_bwt,
			   fasta_ann,
			   fasta_amb,
			   fasta_pac) | set { gridss_output }
		}

		GRIDSS_FILTER_VCF(gridss_output)

		VEP_SV(GRIDSS_FILTER_VCF.out, fasta, fasta_fai)

		if (params.virusbreakend == true) {
			VIRUS_BREAKEND(ch_bam, 
						   fasta,
						   fasta_fai)
		}
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
        read_bam(params.samplesheet) | set { ch_bam }

        if (params.manta == true) {
            manta(ch_bam, 
                        params.fasta)
        }
        
        if (params.delly == true) {
            delly(ch_bam, 
                    params.fasta,
                    params.centromeres,
                    params.delly_map)
        }

        if (params.smoove == true) {
            smoove(ch_bam,
                        params.fasta)
        }

        if (params.gridss == true) {
			gridss(ch_bam,
						params.fasta)
        }

        // survivor(manta.out,
        //             delly.out,
        //             smoove.out)
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
