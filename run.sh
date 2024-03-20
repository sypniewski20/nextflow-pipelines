nextflow run main.nf -profile docker -resume \
                     --mode mapping \
                     --bqsr false \
                     --cohort_mode false \
                     --exome false \
                     --samplesheet /home/mateuszsypniewski/S8436/samplesheet.tsv \
                     --outfolder /data/output/S8436 \
                     --fasta /data/references/fasta/GRCh38 \
                     --singularity /data/references/singularity