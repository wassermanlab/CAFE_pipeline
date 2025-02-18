// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SNV calling, split the vcf by chr
//	split the multiallelic variants as one per line

process SV_split_vcf_by_chr {
        tag "${chr}"

	input :
	file vcf_file
	val assembly
	val batch
	val run
	each chr
	val var_type

	output :
	path '*.vcf.gz', emit : vcf_onechr	

	script :
	"""
	# Unload bcchr, and load cvmfs
        # unload_bcchr
        source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
        # load cvmfs
        source /cvmfs/soft.computecanada.ca/config/profile/bash.sh

        module load StdEnv/2020
        module load vcftools
	module load bcftools

	vcftools --gzvcf ${vcf_file} --chr ${chr} --recode --recode-INFO-all --out ${vcf_file.simpleName}_${chr}
	#Remove duplicated SV in vcf that appears because of paragraph?
	bcftools norm --rm-dup all -O z -o ${vcf_file.simpleName}_${chr}.vcf.gz ${vcf_file.simpleName}_${chr}.recode.vcf
	"""
}
