#R script
#Created by Solenne Correard in December 2021
#Owned by the Silent Genomes Project Activity 3 team
#Developped to build the IBVL, a background variant library

library(vcfR)
library(dplyr)
library(stringr)
library(tidyr)

###################
# Goal : Create the tables as expected in the SQL
###################

#Read the arguments from the nextflow process calling the Rscript named 'MT_heteroplasmy_bin.nf'
args <- commandArgs(trailingOnly = TRUE)

#Define assembly from the args
assembly=(args[1])
var_type=(args[4])

algorithm="MELT"

#Read the vcf file with annotation and variant frequencies
frequ_annot_file=read.vcfR(args[2])
chromosome=frequ_annot_file@fix[1,1]

table_frequ_SV=data.frame()
table_annot_SV=data.frame()
table_sv_consequence=data.frame()

vcf = frequ_annot_file@fix
chr = vcf[,"CHROM"]
pos = as.numeric(vcf[,"POS"])
variant = vcf[,"ID"]
ref = vcf[,"REF"]
alt = vcf[,"ALT"]
quality = vcf[,"QUAL"]

info = vcf[,"INFO"]

#Info from Hail
#Extract table of information
length = sapply(strsplit(info, "SVLEN="),"[[",2)
length = as.data.frame(sapply(strsplit(length, ";"),"[[",1))
show("length ok")


#TSD : Provides the sequence of the Target Site Duplication (TSD) for each insertion. Will be ‘null’ if was not able to be determined.
TSD = sapply(strsplit(info, "TSD="),"[[",2)
TSD =  as.data.frame(sapply(strsplit(TSD, ";"),"[[",1))
show("TSD ok")

# ASSESS : Provides an accuracy assessment, and the information used, to determine the breakpoint of each insertion. Please read our paper to learn more above individual scores (see Citing MELT).
ASSESS = sapply(strsplit(info, "ASSESS="),"[[",2)
ASSESS =  as.data.frame(sapply(strsplit(ASSESS, ";"),"[[",1))
show("ASSESS ok")

#Internal : Provides information on if the insertion is internal to a gene in the provided reference annotation file. The first piece of information is the gene name; the second is the location within that gene (either exon, intron, 5_UTR, or 3_UTR).
INTERNAL = sapply(strsplit(info, "INTERNAL="),"[[",2)
INTERNAL =  as.data.frame(sapply(strsplit(INTERNAL, ";"),"[[",1))
show("INTERNAL ok")

SVTYPE = sapply(strsplit(info, "SVTYPE="),"[[",2)
SVTYPE =  as.data.frame(sapply(strsplit(SVTYPE, ";"),"[[",1))
show("SVTYPE ok")	

#Frequ_table
af_tot_xx_xy = (sapply(strsplit(info, "AF_tot_XX_XY="),"[[",2))
af_tot_xx_xy = (sapply(strsplit(af_tot_xx_xy, ";"),"[[",1))

ac_tot_xx_xy = (sapply(strsplit(info, "AC_tot_XX_XY="),"[[",2))
ac_tot_xx_xy = (sapply(strsplit(ac_tot_xx_xy, ";"),"[[",1))

an_tot_xx_xy = (sapply(strsplit(info, "AN_tot_XX_XY="),"[[",2))
an_tot_xx_xy = (sapply(strsplit(an_tot_xx_xy, ";"),"[[",1))

hom_tot_xx_xy = (sapply(strsplit(info, "hom_tot_XX_XY="),"[[",2))
hom_tot_xx_xy = (sapply(strsplit(hom_tot_xx_xy, ";"),"[[",1))

af_tot = as.data.frame(sapply(strsplit(af_tot_xx_xy, ","),"[[",1))
ac_tot = as.data.frame(sapply(strsplit(ac_tot_xx_xy, ","),"[[",1))
an_tot = as.data.frame(sapply(strsplit(an_tot_xx_xy, ","),"[[",1))
hom_tot = as.data.frame(sapply(strsplit(hom_tot_xx_xy, ","),"[[",1))

af_xx = as.data.frame(sapply(strsplit(af_tot_xx_xy, ","),"[[",2))
ac_xx = as.data.frame(sapply(strsplit(ac_tot_xx_xy, ","),"[[",2))
an_xx = as.data.frame(sapply(strsplit(an_tot_xx_xy, ","),"[[",2))
hom_xx = as.data.frame(sapply(strsplit(hom_tot_xx_xy, ","),"[[",2))

af_xy = as.data.frame(sapply(strsplit(af_tot_xx_xy, ","),"[[",3))
ac_xy = as.data.frame(sapply(strsplit(ac_tot_xx_xy, ","),"[[",3))
an_xy = as.data.frame(sapply(strsplit(an_tot_xx_xy, ","),"[[",3))
hom_xy = as.data.frame(sapply(strsplit(hom_tot_xx_xy, ","),"[[",3))

show("frequ ok")

##Generate table with info and annot
#VEP annotation
vep_annot = as.data.frame(sapply(strsplit(info, "CSQ"),"[[",2))

#Split annot as one varaint can have several annot
#Create large table
all_info = cbind(chr, pos, variant, ref, alt, quality, af_tot, ac_tot, an_tot, hom_tot, af_xx, ac_xx, an_xx, hom_xx, af_xy, ac_xy, an_xy, hom_xy, length, TSD, ASSESS, INTERNAL, SVTYPE, vep_annot)
colnames(all_info) = c ("chr", "pos", "variant", "ref", "alt", "quality", "af_tot", "ac_tot", "an_tot", "hom_tot", "af_xx", "ac_xx", "an_xx", "hom_xx", "af_xy", "ac_xy", "an_xy", "hom_xy", "length", "TSD", "ASSESS", "INTERNAL", "SVTYPE", "vep_annot")

#Creates one row per annotation
all_info_split = separate_rows(all_info, vep_annot, sep = ",")

vep_annot = all_info_split$vep_annot

#Split the VEP annotation (seperated by |)
vep_annot_split = do.call(rbind.data.frame, strsplit(vep_annot, "\\|"))
colnames(vep_annot_split) = c("Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "TSL", "REFSEQ_MATCH", "SOURCE", "REFSEQ_OFFSET", "GIVEN_REF", "USED_REF", "BAM_EDIT", "SIFT", "PolyPhen", "HGVS_OFFSET", "CLIN_SIG", "SOMATIC", "PHENO", "VAR_SYNONYMS", "CADD_PHRED", "CADD_RAW", "SpliceAI_pred_DP_AG", "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")
        #, "SpliceAI_pred_SYMBOL") 

#Create large table
all_info_complete = cbind(all_info_split$chr, all_info_split$pos, all_info_split$variant, all_info_split$ref, all_info_split$alt, all_info_split$quality, all_info_split$af_tot, all_info_split$ac_tot, all_info_split$an_tot, all_info_split$hom_tot, all_info_split$af_xx, all_info_split$ac_xx, all_info_split$an_xx, all_info_split$hom_xx, all_info_split$af_xy, all_info_split$ac_xy, all_info_split$an_xy, all_info_split$hom_xy, all_info_split$length, all_info_split$TSD, all_info_split$ASSESS, all_info_split$INTERNAL, all_info_split$SVTYPE, vep_annot_split$Allele, vep_annot_split$Consequence, vep_annot_split$IMPACT, vep_annot_split$SYMBOL, vep_annot_split$Gene, vep_annot_split$Feature_type, vep_annot_split$Feature, vep_annot_split$BIOTYPE, vep_annot_split$EXON, vep_annot_split$INTRON, vep_annot_split$HGVSc, vep_annot_split$HGVSp, vep_annot_split$cDNA_position, vep_annot_split$CDS_position, vep_annot_split$Protein_position, vep_annot_split$Amino_acids, vep_annot_split$Codons, vep_annot_split$Existing_variation, vep_annot_split$DISTANCE, vep_annot_split$STRAND, vep_annot_split$FLAGS, vep_annot_split$VARIANT_CLASS, vep_annot_split$SYMBOL_SOURCE, vep_annot_split$HGNC_ID, vep_annot_split$TSL, vep_annot_split$REFSEQ_MATCH, vep_annot_split$SOURCE, vep_annot_split$REFSEQ_OFFSET, vep_annot_split$GIVEN_REF, vep_annot_split$USED_REF, vep_annot_split$BAM_EDIT, vep_annot_split$SIFT, vep_annot_split$PolyPhen, vep_annot_split$HGVS_OFFSET, vep_annot_split$CLIN_SIG, vep_annot_split$SOMATIC, vep_annot_split$PHENO, vep_annot_split$VAR_SYNONYMS, vep_annot_split$CADD_PHRED, vep_annot_split$CADD_RAW, vep_annot_split$SpliceAI_pred_DP_AG, vep_annot_split$SpliceAI_pred_DP_AL, vep_annot_split$SpliceAI_pred_DP_DG, vep_annot_split$SpliceAI_pred_DP_DL, vep_annot_split$SpliceAI_pred_DS_AG, vep_annot_split$SpliceAI_pred_DS_AL, vep_annot_split$SpliceAI_pred_DS_DG, vep_annot_split$SpliceAI_pred_DS_DL)	  
			  
colnames(all_info_complete) = c("chr", "pos", "variant", "ref", "alt", "quality", "af_tot", "ac_tot", "an_tot", "hom_tot", "af_xx", "ac_xx", "an_xx", "hom_xx", "af_xy", "ac_xy", "an_xy", "hom_xy", "length", "TSD", "ASSESS", "INTERNAL", "SVTYPE", "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "TSL", "REFSEQ_MATCH", "SOURCE", "REFSEQ_OFFSET", "GIVEN_REF", "USED_REF", "BAM_EDIT", "SIFT", "PolyPhen", "HGVS_OFFSET", "CLIN_SIG", "SOMATIC", "PHENO", "VAR_SYNONYMS", "CADD_PHRED", "CADD_RAW", "SpliceAI_pred_DP_AG", "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")
        #, "SpliceAI_pred_SYMBOL")  

all_info = as.data.frame(all_info_complete)
	
all_info$MEI_type=paste0("MEI,", all_info$SVTYPE)
all_info$algorithm = algorithm

all_info$start = as.numeric(all_info$pos)
all_info$length = as.numeric(all_info$length)
all_info$end = (all_info$start + all_info$length)

show("end")
show(all_info$end)

# variant_ID, type, length, chr, pos, ref, alt, cadd_score, cadd_interpr, dbsnp_id, dbsnp_url, UCSC_url, ensembl_url, clinvar_url, gnomad_url
# UCSC URL : https://genome.ucsc.edu/cgi-bin/hgTracks?db=<assembly>&highlight=<assembly>.chrM%3A<pos>-<pos>&position=chrM%3A<pos-25>-<pos+25> / hg38 or hg19 db=hg38&highlight=hg38.chrM%3A8602-8602&position=chrM%3A8577-8627
# Ensembl_url : https://uswest.ensembl.org/Homo_sapiens/Location/View?r=<chr>%3A<pos-25>-<pos+25> // r=17%3A63992802-64038237
# clinvar URL : https://www.ncbi.nlm.nih.gov/clinvar/variation/<VCV>/ // 692920
# gnomad_URL : https://gnomad.broadinstitute.org/variant/M-<pos>-<ref>-<alt>?dataset=gnomad_r3 // M-8602-T-C
  

# UCSC URL : https://genome.ucsc.edu/cgi-bin/hgTracks?db=<assembly>&highlight=<assembly>.chrM%3A<pos>-<pos>&position=chrM%3A<pos-25>-<pos+25> / hg38 or hg19 db=hg38&highlight=hg38.chrM%3A8602-8602&position=chrM%3A8577-8627
all_info$ucsc_url=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=",assembly, "&highlight=", assembly, ".chr", all_info$chr, "%3A", all_info$start, "-", all_info$end, "&position=chr", all_info$chr, "%3A", all_info$start-(0.25*all_info$length), "-", all_info$end+(0.25*all_info$length))
  
	#Could be added in V2, not part of the current SQL
	# Ensembl_url : https://uswest.ensembl.org/Homo_sapiens/Location/View?r=<chr>%3A<pos-25>-<pos+25> // r=17%3A63992802-64038237
	#if (assembly=="GRCh38") {
	#	ensembl_url=paste0("https://uswest.ensembl.org/Homo_sapiens/Location/View?r=", chr, "%3A", AVG_START-(0.25*AVG_LEN), "-", AVG_END+(0.25*AVG_LEN))
	#} else if (assembly=="GRCh37") {
	#	ensembl_url=paste0("https://grch37.ensembl.org/Homo_sapiens/Location/View?r=", chr, "%3A", AVG_START-(0.25*AVG_LEN), "-", AVG_END+(0.25*AVG_LEN))
	#}  
    
	#Could be added in V2, not part of the current SQL
	# dbsnp_id : From annotation file (SNV_annot_i), "Existing_variation" column
	#if (grepl("rs", SV_annot_i$Existing_variation)) {
	#	dbsnp_id = gsub(",.*$", "", SV_annot_i$Existing_variation)
	#	dbsnp_URL : https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=<rs_number> // rs1556423501
	#	dbsnp_url=paste0("https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=", dbsnp_id)
	#} else {
	#	dbsnp_id="NA"
	#	dbsnp_url="NA"
	#}

#For SV V1, the gnomAD URL will point to the region view in gnomAD SV
#gnomAD SV is currently only avaialble in GRCh37
#https://gnomad.broadinstitute.org/region/22-18738600-18747000?dataset=gnomad_sv_r2_1
all_info$gnomad_id = "TBD"
all_info$gnomad_url_region<- with(all_info, ifelse(assembly=="GRCh37", paste0("https://gnomad.broadinstitute.org/region/", all_info$chr, "-", all_info$AVG_START, "-", all_info$AVG_END, "?dataset=gnomad_sv_r2_1"), "NA"))	

#Replace empty cells by "NA"
all_info[all_info == ""] <- "NA"
  
### Create tables
# SV_IBVL_frequency
# Variant ID, AF_tot, AF_XX, AF_XY, AC_tot, AC_XX, AC_XY, AN_tot, AN_XX, AN_XY, Hom_alt_tot, Hom_alt_XX, Hom_alt_XY, qual
table_frequ_SV = cbind(all_info$variant, all_info$af_tot, all_info$af_xx, all_info$af_xy, all_info$ac_tot, all_info$ac_xx, all_info$ac_xy, all_info$an_tot, all_info$an_xx, all_info$an_xy, all_info$hom_tot, all_info$hom_xx, all_info$hom_xy, all_info$quality)
colnames(table_frequ_SV) = c("variant", "af_tot", "af_xx", "af_xy", "ac_tot", "ac_xx", "ac_xy", "an_tot", "an_xx", "an_xy", "hom_tot", "hom_xx", "hom_xy", "quality")
table_frequ_SV = unique(table_frequ_SV)
write.table(table_frequ_SV, file=paste0("genomic_ibvl_frequencies_", var_type, "_", chromosome, ".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

show("IBVL table ok")

# SV_annotation (svs)
# variant_ID, chr1, chr1_pos1 (start), chr1_po2 (end), type, length, algorithm, ucsc_url, gnomad_id, gnomad_url
table_annot_SV = cbind(all_info$variant, all_info$chr, all_info$start, all_info$end, all_info$MEI_type, all_info$length, all_info$algorithm, all_info$ucsc_url, all_info$gnomad_id, all_info$gnomad_url)
colnames(table_annot_SV) = c("variant", "chr1", "chr1_pos1", "chr1_pos2", "sv_type", "sv_length", "algorithm", "ucsc_url", "gnomad_id", "gnomad_url")
table_annot_SV = unique(table_annot_SV)
write.table(table_annot_SV, file=paste0("svs_", var_type, "_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

show("table annot ok")

#sv_consequences table
#gene, variant, consequence (intronic, intergenic, etc)
#No transcript associated to SV, so no ensembl or Refseq
# If there is several consequences on the same line (separrated by a coma), create one line per consequence
table_sv_consequence = cbind(all_info$SYMBOL, all_info$variant, all_info$Consequence)
colnames(table_sv_consequence) = c("gene", "variant", "consequence")
table_variant_consequence_split = as.data.frame(table_sv_consequence)
table_variant_consequence_split=separate_rows(table_variant_consequence_split, consequence, sep = "&")
table_variant_consequence_split = unique(table_variant_consequence_split)
table_variant_consequence_split = as.data.frame(table_variant_consequence_split)
table_variant_consequence_split_noNA = table_variant_consequence_split[!grepl("NA",table_variant_consequence_split$gene), ]
write.table(table_variant_consequence_split_noNA, file=paste0("sv_consequences_", var_type, "_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

show("sv consequ ok")

#Variants
variants_table=unique(cbind(all_info$variant, "SV"))
colnames(variants_table)=c("variant_id", "var_type")
variants_table = unique(variants_table)
write.table(variants_table, file=paste0("variants_sv_", var_type, "_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

show("variants ok")

# Gene table
# Short name, (NCBI gene Id, may not be necessary anymore)
gene_table=as.data.frame(unique(all_info$SYMBOL))
colnames(gene_table)=c("short_name")
gene_table = as.data.frame(gene_table)
gene_table = unique(gene_table)
gene_table_noNA = gene_table[!grepl("NA",gene_table$short_name), ]
gene_table_noNA = as.data.frame(gene_table_noNA)
colnames(gene_table_noNA)=c("short_name")
write.table(gene_table_noNA, file=paste0("genes_sv_", var_type, "_", chromosome, ".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

show("gene ok")

