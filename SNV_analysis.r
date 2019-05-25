install.packages("TCGAretriever") #CRAN package
BiocManager::install("TCGAbiolinks")#bioconductor packagea
BiocManager::install("TCGAWorkflow")
BiocManager::install("BSgenome")# for sequence extraction
BiocManager::install("maftools")
library(maftools)
library(BiocManager)
install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
getwd()
TCGA.BLCA <- read.delim("C:/Users/user/Desktop/INTBIO_Project/Data/SNV/BC/tcga_BC.maf.gz/TCGA.BLCA.varscan.4a755399-e5b5-4d0e-b1e3-dae24b81590e.DR-10.0.somatic.maf", header=FALSE, comment.char="#")
clinical <- read.delim("~/R/win-library/3.5/maftools/extdata/clinical.tsv")

BLCA.maf = system.file('extdata', 'TCGA.BLCA.varscan.4a755399-e5b5-4d0e-b1e3-dae24b81590e.DR-10.0.somatic.maf', package = 'maftools') 
BLCA.clin = system.file('extdata', 'clinical.tsv', package = 'maftools') 

TCGA.BLCA$Tumor_Sample_Barcode= paste (TCGA.BLCA$V16, TCGA.BLCA$V17)
colnames(TCGA.BLCA)
colnames(clinical)
colnames(clinical)= c("Tumor_Sample_Barcode","submitter_id","project_id","gender","year_of_birth","race","ethnicity","year_of_death","classification_of_tumor","last_known_disease_status","primary_diagnosis","tumor_stage","age_at_diagnosis","vital_status","morphology","days_to_death","days_to_last_known_disease_status","days_to_recurrence","tumor_grade","tissue_or_organ_of_origin","days_to_birth","progression_or_recurrence","prior_malignancy","site_of_resection_or_biopsy","days_to_last_follow_up","therapeutic_agents","treatment_intent_type" )

BLCA = read.maf(maf = BLCA.maf)
getSampleSummary(BLCA)
getGeneSummary(BLCA)
getClinicalData(BLCA)
getFields(BLCA)
#mafsumary plot
plotmafSummary(maf = BLCA, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#top10 mutated genes
oncoplot(maf = BLCA, top = 10, fontSize = 12)
#enrichment of known oncogenic pathways
OncogenicPathways(BLCA)
#de-novo signature analysis (Extract single 5â and 3â bases flanking the mutated site)
library(maftools)
BLCA.tnm= trinucleotideMatrix(BLCA, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19', prefix = NULL, add = TRUE, ignoreChr = NULL, useSyn = TRUE, fn = NULL)
#BSgenome::available.genomes()
#getBSgenome("BSgenome.Hsapiens.UCSC.hg19", masked=FALSE)


extractSignatures(mat, n = NULL, nTry = 6, plotBestFitRes = FALSE,parallel = NULL, pConstant = NULL)#result in sig_res/Decompose a matrix of 96 substitution classes into n signatures.
signatureEnrichment(maf, sig_res, minMut = 5, useCNV = FALSE,fn = NULL)#result list containing p-values/Performs k-means clustering to assign signature to samples and performs enrichment analysis.
                                    

#Exact tests to detect mutually exclusive, co-occuring and altered genesets.
somaticInteractions(maf= BLCA, top = 25, genes = c("ABCG5","ADCY1","BAALC","CDA","CEACAM5","CRNN","CXCL6","FLG2","FOLR1","GJB2","HEPHL1","IL17RB","KRT1","KRTDAP","LOC107984749","LOC729867","MAGEA6","MAGEA8","MTRNR2L1","MUC13","MUC21","MUCL3","MYCN","MYH2","NKX2-1","PRAME","PSORS1C3","RHCG","S100A7","SCEL","SLCO1B3","SPDEF","SULT4A1","TCN1","TGM5","TM4SF19","WNK4"), pvalue = c(0.05, 0.01), returnAll = FALSE, findPathways = TRUE, kMax = 37,fontSize = 0.8, verbose = TRUE)

##===================
all= TCGA.BLCA
candidateGenes= read.table("C:/Users/user/Desktop/INTBIO_Project/BLadder_Results/Normalvsall - Copy.txt", quote="\"", comment.char="")
subset=candidateGenes  
x= all[all$V1 %in% subset$V1, ]
colnames(x)= c("Hugo_Symbol","Entrez_Gene_Id","Center","NCBI_Build","Chromosome","Start_Position","End_Position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","dbSNP_RS","dbSNP_Val_Status","Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode","Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2","Tumor_Validation_Allele1","Tumor_Validation_Allele2","Match_Norm_Validation_Allele1","Match_Norm_Validation_Allele2", "Verification_Status","Validation_Status","Mutation_Status","Sequencing_Phase","Sequence_Source","Validation_Method","Score","BAM_File","Sequencer","Tumor_Sample_UUID","Matched_Norm_Sample_UUID","HGVSc,HGVSp,HGVSp_Short","Transcript_ID","Exon_Number","t_depth","t_ref_count","t_alt_count","n_depth","n_ref_count","n_alt_count","all_effects","Allele","Gene","Feature","Feature_type","One_Consequence","Consequence","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","ALLELE_NUM","DISTANCE","TRANSCRIPT_STRAND","SYMBOL","SYMBOL_SOURCE","HGNC_ID","BIOTYPE","CANONICAL","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","RefSeq","SIFT","PolyPhen","EXON","INTRON","DOMAINS","GMAF","AFR_MAF","AMR_MAF","ASN_MAF","EAS_MAF","EUR_MAF","SAS_MAF","AA_MAF","EA_MAF","CLIN_SIG","SOMATIC","PUBMED","MOTIF_NAME","MOTIF_POS","HIGH_INF_POS","MOTIF_SCORE_CHANGE","IMPACT","PICK","VARIANT_CLASS","TSL","HGVS_OFFSET","PHENO", "MINIMISED","ExAC_AF","ExAC_AF_Adj","ExAC_AF_AFR","ExAC_AF_AMR","ExAC_AF_EAS","ExAC_AF_FIN","ExAC_AF_NFE","ExAC_AF_OTH","ExAC_AF_SAS", "GENE_PHENO","FILTER","CONTEXT","src_vcf_id","tumor_bam_uuid","normal_bam_uuid","case_id","GDC_FILTER","COSMIC,MC3_Overlap","GDC_Validation_Status")
myvars <- c("Hugo_Symbol", "dbSNP_RS")
newdata <- x[myvars]


#from dgidb database
Drug_interaction= read.delim("C:/Users/user/Desktop/INTBIO_Project/BLadder_Results/dgidb_export_2019-03-04.tsv")
data=write.xlsx(x, "c:/mydata.xlsx")
