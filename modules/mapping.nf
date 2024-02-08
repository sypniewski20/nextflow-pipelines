// process BWA_MAP_READS {
// 	tag "${sample}"
// 	label 'gatk'
// 	label 'mem_64GB'
// 	label 'core_36'
// 	input:
// 		tuple val(sample), path(read_1), path(read_2)
//         path(fasta)
// 	output:
// 		tuple val(sample), path("${sample}_pre.bam"), path("${sample}_pre.bam.bai")
// 	script:
// 		"""

// 		READ_GROUP_NAME=\$( zgrep @ ${read_1} | \
// 					head -n 1 | \
// 					cut -d :  -f1 | \
// 					sed 's/@//g' )

// 		bwa mem -t ${task.cpus} \
// 			${fasta}/${fasta}.fa \
// 			${read_1} \
// 			${read_2} | \
// 		samtools view --reference ${fasta}/${fasta}.fa \
// 					  --threads ${task.cpus} \
// 					  --bam | \
// 		samtools sort -@ ${task.cpus} \
// 					  -O bam \
// 					  -o /dev/stdout | \
// 		gatk AddOrReplaceReadGroups \
// 					  -I /dev/stdin \
// 					  -O  ${sample}_pre.bam \
// 					  -RGID \$READ_GROUP_NAME \
// 					  -RGLB lib1 \
// 					  -RGPL ILLUMINA \
// 					  -RGPU unit1 \
// 					  -RGSM ${sample} \
// 					  --CREATE_INDEX true

// 		samtools index -@ ${task.cpus} ${sample}_pre.bam

// 		gatk ValidateSamFile -I  ${sample}_pre.bam -R ${fasta}/${fasta}.fa


// 		"""
// }

process BWA_MAP_READS {
	tag "${sample}"
	label 'gatk'
	label 'mem_64GB'
	label 'core_36'
	input:
		tuple val(sample), path(read_1), path(read_2)
        path(fasta)
	output:
		tuple val(sample), path("${sample}_pre.bam"), path("${sample}_pre.bam.bai")
	script:
		"""

		READ_GROUP_NAME=\$( zgrep @ ${read_1} | \
					head -n 1 | \
					cut -d :  -f1 | \
					sed 's/@//g' )

		tags=\$(echo "@RG\\tID:\${READ_GROUP_NAME}\\tSM:${sample}\\tPL:ILLUMINA\\tLB:lb1")

		bwa mem -t ${task.cpus} \
			-R \${tags} \
			${fasta}/${fasta}.fa \
			${read_1} \
			${read_2} | \
		samtools view --reference ${fasta}/${fasta}.fa \
					  --threads ${task.cpus} \
					  --bam | \
		samtools sort -@ ${task.cpus} \
					  -O bam \
					  -o ${sample}_pre.bam

		samtools index -@ ${task.cpus} ${sample}_pre.bam

		gatk ValidateSamFile -I  ${sample}_pre.bam \
							 -MODE	SUMMARY


		"""
}


process BASE_RECALIBRATOR {
	publishDir "${params.outfolder}/${params.runID}/BAM", mode: 'symlink', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		path(fasta)
		path(interval_list)
		path(snv_resource)
		path(cnv_resource)
		path(sv_resource)
	output:
		tuple val(sample), path("${bam}"), path("${sample}_recal_data.table")
	script:
		"""

		for res in ${snv_resource} ${cnv_resource} ${sv_resource}; do
			
			ln -sf \$( echo "\$( realpath \${res} ).tbi" ) .
			
		done

		gatk BaseRecalibrator \
		-I ${bam} \
		-R ${fasta}/${fasta}.fa \
		-L ${interval_list} \
		--known-sites ${snv_resource} \
		--known-sites ${cnv_resource} \
		--known-sites ${sv_resource} \
		-O ${sample}_recal_data.table
		"""
}

process APPLY_BQSR {
	publishDir "${params.outfolder}/${params.runID}/BAM", mode: 'symlink', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_8'
	input:
		tuple val(sample), path(bam), path(recal_table)
		path(interval_list)
		path(fasta)
	output:
		tuple val(sample), path("${sample}_recal.bam"), path("${sample}_recal.bam.bai")
	script:
		"""
		gatk ApplyBQSR \
		-R ${fasta}/${fasta}.fa \
		-L ${interval_list} \
		-I ${bam} \
		--bqsr-recal-file ${recal_table} \
		-O ${sample}_recal.bam

		samtools index -@ ${task.cpus} ${sample}_recal.bam

		gatk ValidateSamFile -I ${sample}_recal.bam \
							 -M SUMMARY

		"""
}

process GATK_MARK_DUPLICATES {
	publishDir "${params.outfolder}/${params.runID}/BAMQC", pattern: "${sample}.dup_metrics.*", mode: 'copy', overwrite: true
	publishDir "${params.outfolder}/${params.runID}/BAM", pattern: "${sample}_recal_dupfiltered.bam*", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_1'
	input:
		tuple val(sample), path(bam), path(bai)
		path(interval_list)
	output:
		tuple val(sample),  path("${sample}_recal_dupfiltered.bam"), path("${sample}_recal_dupfiltered.bam.bai"), emit: ch_bam
		path("${sample}.txt"), emit: sample_checkpoint
	script:
		"""
       	gatk MarkDuplicates \
            -I ${bam} \
            -O ${sample}_recal_dupfiltered.bam \
            -M ${sample}.dup_metrics.txt \
			--MAX_RECORDS_IN_RAM 5000000 \
			--CREATE_INDEX true
		
		if \$( samtools flagstat ${sample}_recal_dupfiltered.bam | grep -q "^0 + 0 in total (QC-passed reads + QC-failed reads)\$" ); then
        	echo "${sample}_recal_dupfiltered.bam is empty"
			exit 1
		fi

		echo -e "${sample}\\t${is_germline}\\t${params.outfolder}/${params.runID}/BAM/${sample}_recal_dupfiltered.bam\\t${params.outfolder}/${params.runID}/BAM/${sample}_recal_dupfiltered.bam.bai" > ${sample}.txt
		"""
}

process SAMBAMBA_MARK_DUPLICATES {
	publishDir "${params.outfolder}/${params.runID}/BAMQC", pattern: "${sample}.dup_metrics.*", mode: 'copy', overwrite: true
	publishDir "${params.outfolder}/${params.runID}/BAM", pattern: "${sample}_recal_dupfiltered.*", mode: 'copy', overwrite: true
	tag "${sample}"
	label 'gatk'
	label 'mem_16GB'
	label 'core_8'
	input:
		tuple val(sample), path(bam), path(bai)
	output:
		tuple val(sample),  path("${sample}_recal_dupfiltered.bam"), path("${sample}_recal_dupfiltered.bam.bai"), emit: ch_bam
		path("${sample}.txt"), emit: sample_checkpoint
	script:
		"""

		sambamba markdup \
			--nthreads ${task.cpus} \
			${bam} \
			${sample}_recal_dupfiltered.bam
		
		sambamba index --nthreads ${task.cpus} ${sample}_recal_dupfiltered.bam

		if \$( sambamba flagstat --nthreads ${task.cpus} ${sample}_recal_dupfiltered.bam | grep -q "^0 + 0 in total (QC-passed reads + QC-failed reads)\$" ); then
        	echo "${sample}_recal_dupfiltered.bam is empty"
			exit 1
		fi

		gatk ValidateSamFile -I ${sample}_recal_dupfiltered.bam \
							 -M SUMMARY

		echo -e "${sample}\\t${params.outfolder}/${params.runID}/BAM/${sample}_recal_dupfiltered.bam\\t${params.outfolder}/${params.runID}/BAM/${sample}_recal_dupfiltered.bam.bai" > ${sample}.txt
		
		"""
}

