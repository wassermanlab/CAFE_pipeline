// Nextflow process
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics:
// Fix the bam files adding tags to allow faster processing by some future steps


process samtools_fixmate {
	label 'conda_annotate'
	tag "${bam.SimpleName}"

	input :
	file bam
	file bai
	val assembly
	val batch
	val run

	output :
	path '*_fixmate_ordered.bam', emit: samples_fixmate_bam
	path '*_fixmate_ordered.bam.bai', emit: samples_fixmate_bam_index

	script:
	"""
	sample_name=\$(echo ${bam.simpleName} | cut -d _ -f 1)
		# Resort the bam file by query name for samtools fixmate (coordiante-sorted bam does not work)
		samtools sort -n -O BAM -@ 20  ${bam} > ${bam.SimpleName}_nsorted.bam

		# Samtools fixmate will add MQ tags
		samtools fixmate -m -O BAM -@ 20  ${bam.SimpleName}_nsorted.bam  ${bam.SimpleName}_fixmate.bam

		# Now sort bam file by coordinates to resume the pipeline 
		samtools sort  -@ 20  ${bam.SimpleName}_fixmate.bam -o ${bam.SimpleName}_fixmate_ordered.bam

		# index the sorted bam file
		samtools index  -@ 20  ${bam.SimpleName}_fixmate_ordered.bam
		rm *_nsorted.bam* *_fixmate.bam
	
	"""
}

