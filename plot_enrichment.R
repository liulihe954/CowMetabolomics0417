setwd("/Users/liulihe95/Desktop/GuanYang")
#BiocManager::install("clusterProfiler")
library(clusterProfiler);library(ggplot2);library(dplyr);library(magrittr)
require("biomaRt");require("gage")
#use mart and select dataset: here use names as connection
mart <- useMart("ENSEMBL_MART_MOUSE")
mart <- useDataset("mmusculus_gene_ensembll", mart)
mart <- useDataset("m129s1svimj_gene_ensembl", mart)

#listDatasets(mart) # show all dataset available
test1 = read.csv("sig1.csv",header = T,sep = ",")
test2 = read.csv("sig2.csv",header = T,sep = ",")
signif_CD4T = test1[!is.na(test1$Symbol),3]
signif_CD8T = test2[!is.na(test2$Symbol),3]
# convert names to ENTERid (either way) to get annotation dataset
annot_CD4T = getBM(attributes=c("ensembl_gene_id", "entrezgene","external_gene_name"),
                    filters="external_gene_name",
                     values=signif_CD4T,
                     mart=mart)
annot_CD8T = getBM(attributes=c("ensembl_gene_id", "entrezgene","external_gene_name"),
                   filters="external_gene_name",
                   values=signif_CD8T,
                   mart=mart)
# extract enterID from the annotation
enterID_CD4T = annot_CD4T$entrezgene;enterID_CD8T = annot_CD8T$entrezgene
# get enrichment results
enrich_CD4T <- enrichKEGG(gene = enterID_CD4T, organism = 'mmu', qvalueCutoff = 0.05, pvalueCutoff = 0.05)
enrich_CD8T <- enrichKEGG(gene = enterID_CD8T, organism = 'mmu', qvalueCutoff = 0.05, pvalueCutoff = 0.05)

### need massage here
enrich_CD4T$GeneRatio;enrich_CD4T$BgRatio
overlap1 = sub('/.*', '',enrich_CD4T$GeneRatio)
total1 = sub('/.*', '',enrich_CD4T$BgRatio)

enrich_CD4T_add = cbind(data.frame(enrich_CD4T),total1 = as.numeric(total1),overlap1 = as.numeric(overlap1))
test_enrich_CD4T = enrich_CD4T$GeneRatio
enrich_CD4T_add 

pdf("test.guan1.pdf",paper = "A4")
enrich_CD4T_add %>%
  top_n(7, wt= -pvalue) %>%
  mutate(hitsPerc=overlap1*100/total1) %>% ## signi genes, v1 = all genes in the go.
  ggplot(aes(x=hitsPerc,
             y=Description,
             colour=pvalue,
             size=overlap1)) +
  xlim(0,22)+
  geom_point() +
  theme_gray()+
  theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  #expand_limits(x=0) +
  labs(x="Hits (%)", y="Kegg term CD4T", colour="p value", size="Count")
dev.off()




lidat = read.csv("lidat.txt",header = T,sep = "")

s
substr(enrich_CD8T$GeneRatio,1,2)

#  extract the number of genes
#in bash
#cat li.txtx |awk -F'"' '{print $10}'| egrep -v p | awk -F'/'  '{print $1} ' > total.gens

pdf("kegg.li.pdf")
lp %>%
  top_n(10, wt= -pvalue) %>%
  mutate(hitsPerc= ligenes*100/V1) %>% ## signi genes, v1 = all genes in the go.
  ggplot(aes(x=hitsPerc,
             y= Description,
             colour=pvalue,
             size= ligenes)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term BP", colour="p value", size="Count")



##########################################################################################################################################
## Prepare data for Gene Set Analysis

total.genes = c() # total genes in your dataset
sig.genes = c() # total genes in the non-preserved module

## Analysis bosTau annotation: GO
library(biomaRt)

database = useMart("ensembl")
genome = useDataset("btaurus_gene_ensembl", mart = database)
gene = getBM(c("ensembl_gene_id","external_gene_name","go_id","name_1006"), mart = genome)

dim(gene); length(unique(gene$ensembl_gene_id)); length(unique(gene$go_id))

goName = unique(gene[,c(3,4)])
goName = goName[order(goName$go_id),]
goName = goName[-1,]

GO = goName$go_id
Name = goName$name_1006
genesGO = unique(subset(gene,go_id != "")$ensembl_gene_id)
N = length(total.genes[total.genes %in% genesGO])
S = length(sig.genes[sig.genes %in% genesGO])
out = data.frame(GO=character(),Name=character(),totalG=numeric(),sigG=numeric(),Pvalue=numeric())

for(i in 1:length(GO)){
  gENEs = subset(gene, go_id == GO[i])$ensembl_gene_id 
  m = length(total.genes[total.genes %in% gENEs])
  s = length(sig.genes[sig.genes %in% gENEs])
  M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
  Pval = round(fisher.test(M, alternative ="g")$p.value, digits = 3)
  tmp = data.frame(GO = GO[i], Name = Name[i], totalG = m, sigG = s, Pvalue = Pval)
  out = rbind(out,tmp)}

ot = subset(out,totalG > 4 & Pvalue < 0.05)
final = ot[order(ot$Pvalue),]
colnames(final) = c("GOID","GO Name", "Total Genes", "Significant Genes", "P-value")
