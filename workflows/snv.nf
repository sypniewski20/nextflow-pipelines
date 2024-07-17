nextflow.enable.dsl = 2

println """\
=================================================================================================================
										Run details: 
=================================================================================================================
Run ID              :${params.runID}
Output				:${params.outfolder}/${params.runID}
Call MT				:${params.mity}
=================================================================================================================	
										References
=================================================================================================================
Fasta 				:${params.fasta}
=================================================================================================================
"""
.stripIndent()

include { Read_bam_checkpoint } from '../modules/functions.nf'
include { DEEP_VARIANT_WGS; DEEP_VARIANT_WES; FILTER_SNVS } from "../modules/deepvariant.nf"
include { GET_MT; MT_CALL; MT_REPORT; MT_MERGE} from "../modules/mity.nf"
include { VEP_SNV } from "../modules/annotations.nf"


workflow read_bam {
	take:
		ch_checkpoint
	main:
		Read_bam_checkpoint(ch_checkpoint)
	emit:
		ch_bam
}

workflow deepvariant {
	take:
		ch_bam
		fasta
	main:
		fasta_fai = file(fasta+'.fai')
		if( params.exome == false ) {
			DEEP_VARIANT_WGS(ch_bam, 
							fasta,
							fasta_fai)
			FILTER_SNVS(DEEP_VARIANT_WGS.out, fasta)
		} else {
			DEEP_VARIANT_WES(ch_bam, 
							fasta,
							fasta_fai)
			FILTER_SNVS(DEEP_VARIANT_WES.out, fasta)
		}

		VEP_SNV(FILTER_SNVS.out, fasta, fasta_fai)

		if (params.mity == true) {
			MT_CALL(ch_bam)
			MT_REPORT(MT_CALL.out)
			MT_MERGE(MT_CALL.out, FILTER_SNVS.out)

			VEP_SNV(MT_MERGE.out, fasta, fasta_fai)

		}

}

//Main workflow

workflow {
	main:
		read_bam(params.samplesheet) | set { ch_bam }

		deepvariant(read_bam.out.ch_bam, 
					params.fasta)
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
