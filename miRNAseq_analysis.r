library(readr)
library(org.Hs.eg.db)
library(DESeq2)
data.path="/Users/Roza/Documents/integrative biology/data/miRNA-seq"
files <- list.files(path = data.path, recursive = T, pattern = "txt")


file = files[1]
file.id = strsplit(file,"/")[[1]][1]
temp <- read.table(file.path(data.path, files[1]), header = T)
mirna.exp = temp
temp = temp[,1:3]
temp = temp[-2]
x=duplicated(temp$miRNA_ID)  
sum(x)
temp.data=temp[2]
temp.data=apply(temp.data,2, as.numeric)
temp.data.agg= aggregate(temp.data, list(temp$miRNA_ID),FUN=sum)
rownames(temp.data.agg) = temp.data.agg[,1]
temp.data.agg = temp.data.agg[-1]
colnames(temp.data.agg)=c(file.id)
mirna.exp = temp.data.agg

for(i in 2: length(files))
{
  
  file = files[i]
  file.id = strsplit(file,"/")[[1]][1]
  temp <- read.table(file.path(data.path, files[i]), header = T)
  temp = temp[,1:3]
  temp = temp[-2]
  x=duplicated(temp$miRNA_ID)  
  sum(x)
  temp.data=temp[2]
  temp.data=apply(temp.data,2, as.numeric)
  temp.data.agg= aggregate(temp.data, list(temp$miRNA_ID),FUN=sum)
  rownames(temp.data.agg) = temp.data.agg[,1]
  temp.data.agg = temp.data.agg[-1]
  colnames(temp.data.agg)=c(file.id)
  mirna.exp = merge(data.frame(temp.data.agg), data.frame(mirna.exp), by = 0, all = TRUE)
  #mirna.exp <- merge(temp.data.agg, mirna.exp, by=0, all=T)
  rownames(mirna.exp) <- mirna.exp$Row.names; mirna.exp$Row.names <- NULL
}

x=duplicated(mirna.exp$row.name)  
sum(x)

pheno <- read_delim("~/Documents/integrative\ biology/data/miRNA-seq/gdc_sample_sheet.2018-12-19.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
table(pheno$`Sample Type`)
pheno.names=names(pheno)
names(pheno)= as.character( sapply ( pheno.names, function(x) gsub(" ",".",x)))
table(pheno$Sample.Type)
file.ids.pheno=pheno$File.ID
index.files=match(file.id,file.ids.pheno)
names(mirna.exp)=pheno$Sample.ID[index.files]


all.normal.samples= pheno[ pheno$Sample.Type %in% c("Solid Tissue Normal"),]$Sample.ID
#normal.samples=all.normal.samples[1:sample.no]

all.tumor.samples= pheno[ pheno$Sample.Type %in% c("Primary Tumor"),]$Sample.ID
#tumor.samples=all.tumor.samples[1:sample.no]

normal.exp=mirna.exp[, names(mirna.exp)%in% all.normal.samples]
tumor.exp=mirna.exp[, names(mirna.exp)%in% all.tumor.samples]

pheno.sub=pheno[pheno$Sample.ID %in% c(all.normal.samples,all.tumor.samples), c("Sample.ID", "Sample.Type")]

exp.sub=cbind(normal.exp,tumor.exp)
exp.sub=apply (exp.sub, 2,as.integer)
rownames(exp.sub)=rownames(normal.exp)


save(mirna.exp,pheno, exp.sub,pheno.sub ,file="miRNA-seq.RDATA")

cond1="Solid Tissue Normal" 
cond2="Primary Tumor"

dds = DESeqDataSetFromMatrix( countData = exp.sub , colData = pheno.sub , design = ~ Sample.Type)
dds.run = DESeq(dds)


