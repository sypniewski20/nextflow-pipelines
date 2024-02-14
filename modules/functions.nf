def Read_samplesheet(samplesheet) {
	Channel
	    .fromPath (samplesheet)
	    .ifEmpty ( "Samplesheet empty." )
	    .splitCsv ( header:true, sep:'\t' )
	    .map { row -> [ row.sampleID, file(row.R1, checkIfExists: true), file(row.R2, checkIfExists: true) ]}
	    .set { ch_samples_initial }
}

def Read_bam_checkpoint(bam_checkpoint_sheet) {
	Channel
		.fromPath (bam_checkpoint_sheet)
	    .ifEmpty ( "Samplesheet empty." )
	    .splitCsv ( header:true, sep:'\t' )
	    .map { row -> [ row.sampleID, file(row.bam, checkIfExists: true), file(row.bai, checkIfExists: true) ]}
	    .set { ch_samples_checkpoint }
}

