setwd("/Users/liulihe95/Desktop/metabolomics_0417/")
library(biomaRt)
library(pathview)
library(KEGGgraph)
library(gage)
##prepare pathway - - - bos taurus
sdb = kegg.gsets(species = "bta", id.type = "kegg", check.new=FALSE)
kegg.gs = sdb$kg.sets[sdb$sigmet.id]

length(sdb$kg.sets)
dim(kegg.gs)


###################################
### dataset prepare for matching 
# up_LF_LCL_ampulla.txt
# up_LS_LCL_isthmus.txt
# up_SF_SCL_ampulla.txt
# up_SF_SCL_isthmus.txt
up_LF_LCL_ampulla = read.csv("up_LF_LCL_ampulla.txt",sep = "\t")
names(up_LF_LCL_ampulla) = c("ensembl_gene_id","Gene Symbol","Base Mean","log2 Fold Change","P value","P adj")
up_SF_SCL_ampulla = read.csv("up_SF_SCL_ampulla.txt",sep = "\t")
names(up_SF_SCL_ampulla) = c("ensembl_gene_id","Gene Symbol","Base Mean","log2 Fold Change","P value","P adj")
up_LF_LCL_isthmus = read.csv("up_LS_LCL_isthmus.txt",sep = "\t")
names(up_LF_LCL_isthmus) = c("ensembl_gene_id","Gene Symbol","Base Mean","log2 Fold Change","P value","P adj")
up_SF_SCL_isthmus = read.csv("up_SF_SCL_isthmus.txt",sep = "\t")
names(up_SF_SCL_isthmus) = c("ensembl_gene_id","Gene Symbol","Base Mean","log2 Fold Change","P value","P adj")

dim(up_LF_LCL_ampulla)
dim(up_SF_SCL_ampulla)
dim(up_LF_LCL_isthmus)
dim(up_SF_SCL_isthmus)



# get mart (database)
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("btaurus_gene_ensembl", mart)

#str(up_SF_SCL_isthmus)

### up_LF_LCL_ampulla.txt
IDs1 = as.vector(up_LF_LCL_ampulla$ensembl_gene_id) #ID in our dataset

annot1 = getBM(attributes=c("ensembl_gene_id","entrezgene_id","external_gene_name"),
               filters="ensembl_gene_id",
               values=IDs1,
               mart=mart)

listAttributes(mart)[,1][grepl("entrezgene",listAttributes(mart)[,1])]

# Problem? Question? Confusion!
length(IDs1)
str(annot1)
# length(!is.na(annot1$external_gene_name))
length(!is.na(annot1$ensembl_gene_id))
length(unique(annot1$ensembl_gene_id))


#annot1 <- annot1[!duplicated(annot1$ensembl_gene_id),]
#length(!is.na(annot1$external_gene_name))

annot1$ensembl_gene_id[duplicated(annot1$ensembl_gene_id)]
annot1[annot1$ensembl_gene_id== "ENSBTAG00000027854",]

### up_LF_LCL_isthmus.txt
IDs2 = as.vector(up_LF_LCL_isthmus$ensembl_gene_id)
annot2 = getBM(attributes=c("ensembl_gene_id", "entrezgene",
                            "external_gene_name"),
               filters="ensembl_gene_id", values=IDs2, mart=mart)
### up_SF_SCL_ampulla.txt
IDs3 = as.vector(up_SF_SCL_ampulla$ensembl_gene_id)
annot3 = getBM(attributes=c("ensembl_gene_id", "entrezgene",
                            "external_gene_name"),
               filters="ensembl_gene_id", values=IDs3, mart=mart)
### up_SF_SCL_isthmus.txt
IDs4 = as.vector(up_SF_SCL_isthmus$ensembl_gene_id)
annot4 = getBM(attributes=c("ensembl_gene_id", "entrezgene",
                            "external_gene_name"),
               filters="ensembl_gene_id", values=IDs4, mart=mart)



###############################################################################
###                 select target pathway and find overlap        #############
##############################################################################
#     1. 
#    "3-Hydroxy Palmitoylcarnitine" 
#       KEYWORD: "Palmitoylcarnitine"
#
#     COMPOUND: C02990 - L-Palmitoylcarnitine
#              map00071  	Fatty acid degradation
#              map01212  	Fatty acid metabolism
#
#     2. 
#     "Phosphatidylcholine with diacyl residue sum C36:6"
#           ("Phosphatidylcholine with acyl-alkyl residue sum C36:4")
#       KEYWORD: "Phosphatidylcholine"
#
#     COMPOUND: C00157  Phosphatidylcholine; Lecithin; Phosphatidyl-N-trimethylethanolamine; 
#                       1,2-Diacyl-sn-glycero-3-phosphocholine; Choline phosphatide;
#                       3-sn-Phosphatidylcholine
#
#             map00564  	Glycerophospholipid metabolism
#             map00590  	Arachidonic acid metabolism
#             map00591  	Linoleic acid metabolism
#             map00592  	alpha-Linolenic acid metabolism
#       x    ? map01100  	Metabolic pathways (not found; maybe too comprehensive?)
#       x  ?  map01110  	Biosynthesis of secondary metabolites (not found; maybe too comprehensive?)
#       x      map04723  	Retrograde endocannabinoid signaling
#       x      map05231  	Choline metabolism in cancer
# ###  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#     COMPOUND: C03631 Oleoylphosphatidylcholine
#              NA
#     COMPOUND: C03873 Linoleoylphosphatidylcholine; 2-Acyl-1-linoleoyl-sn-glycero-3-phosphocholine
#              NA
#     COMPOUND: C03889 Palmitoylphosphatidylcholine; 1-Palmitoyl-2-acyl-sn-glycero-3-phosphocholine
#              NA
#     COMPOUND: C04230 1-Acyl-sn-glycero-3-phosphocholine; 1-Acyl-sn-glycerol-3-phosphocholine;
#                   alpha-Acylglycerophosphocholine; 2-Lysolecithin; 2-Lysophosphatidylcholine;1-Acylglycerophosphocholine
#
#             map00564  	Glycerophospholipid metabolism
#             map05231  	Choline metabolism in cancer
#     
#     COMPOUND: C04233 2-Acyl-sn-glycero-3-phosphocholine; 2-Acylglycero-3-phosphocholine; 1-Lysophosphatidylcholine;
#                    1-Lysolecithin;3-Lysolecithin
#             map00564  	Glycerophospholipid metabolism
#
#           8 unique ones
########################################################################################

all_IDs = list(annot1,annot2,annot3,annot4)
all_dataset = c("up_LF_LCL_ampulla","up_LF_LCL_isthmus","up_SF_SCL_ampulla","up_SF_SCL_isthmus")
all_pathway = c("bta00071 Fatty acid degradation","bta01212 Fatty acid metabolism",
                "bta00564 Glycerophospholipid metabolism", "bta00591 Linoleic acid metabolism",
                "bta00590 Arachidonic acid metabolism","bta00592 alpha-Linolenic acid metabolism",
                "bta00600 Sphingolipid metabolism", "bta04071 Sphingolipid signaling pathway",
                "bta04217 Necroptosis")
selected_pathway = sdb$kg.sets[intersect(all_pathway, names(sdb$kg.sets))]


test_color = unlist(selected_pathway[2])
names(test_color) = NULL
test_color = noquote(test_color)

write.csv(test_color,"test_bta.txt",row.names = F)

overlap = list();names = character();s = 0
for (i in c(1:length(all_IDs))){
  for ( j in c(1:length(all_pathway))){
    s = s + 1
    annot = data.frame(all_IDs[i])
    single_overlap = annot[annot$entrezgene %in% unlist(selected_pathway[j]),]
    overlap[s] = list(single_overlap)
    names(overlap)[s] = paste(all_dataset[i],"---",all_pathway[j])
  }
}
overlap


#each pathway
#    A  / I 
#  L
#  S

# all_dataset = c("up_LF_LCL_ampulla","up_LF_LCL_isthmus","up_SF_SCL_ampulla","up_SF_SCL_isthmus")
Contigency_table_allpath = list()
Overlap_genes = list()
for (i in c(1:length(all_pathway))){
    # i = 1
    #annot = data.frame(all_IDs[i])
    single_overlap1 = annot1[annot1$entrezgene %in% unlist(selected_pathway[i]),]
    single_overlap2 = annot2[annot2$entrezgene %in% unlist(selected_pathway[i]),]
    single_overlap3 = annot3[annot3$entrezgene %in% unlist(selected_pathway[i]),]
    single_overlap4 = annot4[annot4$entrezgene %in% unlist(selected_pathway[i]),]
    
    count1 = nrow(single_overlap1)
    count2 = nrow(single_overlap2)
    count3 = nrow(single_overlap3)
    count4 = nrow(single_overlap4)
    
    overlap1 = intersect(single_overlap1$ensembl_gene_id,single_overlap2$ensembl_gene_id)
    overlap2 = intersect(single_overlap3$ensembl_gene_id,single_overlap4$ensembl_gene_id)
    overlap3 = intersect(single_overlap1$ensembl_gene_id,single_overlap3$ensembl_gene_id)
    overlap4 = intersect(single_overlap2$ensembl_gene_id,single_overlap4$ensembl_gene_id)
    overlap5 = intersect(single_overlap1$ensembl_gene_id,single_overlap4$ensembl_gene_id)
    overlap6 = intersect(single_overlap2$ensembl_gene_id,single_overlap3$ensembl_gene_id)
    
    #overlap_results = paste()
    
    # find overall --- overlap
    #length(overall_lap)
   # overall_lap = c(single_overlap1$ensembl_gene_id,single_overlap3$ensembl_gene_id)[which(c(single_overlap1$ensembl_gene_id,single_overlap3$ensembl_gene_id) == 
                         # c(single_overlap2$ensembl_gene_id,single_overlap4$ensembl_gene_id))]
    #overall_lap = paste()
   # paste(length(overlap5),"/",length(overlap6))
    # test = cat("hello\nworld\n")
    
    contig = matrix(c(count1,count2,length(overlap1),count3,count4,length(overlap2),length(overlap3),length(overlap4),paste(length(overlap5),"/",length(overlap6)))
                    ,byrow = T,3,3)
    rownames(contig) = c("L","S","overlap");colnames(contig) = c("ampulla","isthmus","overlap")
    Contigency_table_allpath[[i]] = contig
    names(Contigency_table_allpath)[i] = paste(names(selected_pathway[i]))
    index = c(length(overlap1),length(overlap2),length(overlap3),length(overlap4))
    combind_overlap = list(overlap1,overlap2,overlap3,overlap4)
    names(combind_overlap) = c(paste(colnames(contig)[1],"---",paste(rownames(contig)[2],"and",rownames(contig)[1])),
                               paste(colnames(contig)[2],"---",paste(rownames(contig)[2],"and",rownames(contig)[1])),
                               paste(rownames(contig)[1],"---",paste(colnames(contig)[1],"and",colnames(contig)[2])),
                               paste(rownames(contig)[2],"---",paste(colnames(contig)[1],"and",colnames(contig)[2])))
    Overlap_genes[[i]] = combind_overlap[which(index >0)]
    names(Overlap_genes)[i] = paste(names(selected_pathway[i]))
}

noquote(Contigency_table_allpath)


which(c(1:3,8) == c(1:4))

## testing
# annot1[annot1$entrezgene %in% unlist(selected_pathway[]),]
# annot2[annot2$entrezgene %in% unlist(selected_pathway[1]),]
# annot3[annot3$entrezgene %in% unlist(selected_pathway[1]),]
# annot4[annot4$entrezgene %in% unlist(selected_pathway[1]),]
# overlap2 = length(intersect(single_overlap3$ensembl_gene_id,single_overlap4$ensembl_gene_id))

## plotting A - up and down (L and S) /// I  up and down (L and S) 
#all_dataset = c("up_LF_LCL_ampulla","up_LF_LCL_isthmus","up_SF_SCL_ampulla","up_SF_SCL_isthmus")
#all_IDs = list(  annot1,              annot2,            annot3,              annot4)
####################################
# Pthway 1-9 
for (i in c(1:length(selected_pathway))){
  #i = 2
  setwd("/Users/liulihe95/Desktop/metabolomics_0417/")
  dir.create("path_ampulla")
  setwd("/Users/liulihe95/Desktop/metabolomics_0417/path_ampulla")
  # ampulla
  pathway_code = substr(names(selected_pathway[i]),4,8)
  index_a = rbind(cbind(annot1[annot1$entrezgene %in% unlist(selected_pathway[i]),],sig = rep(1,nrow(annot1[annot1$entrezgene %in% unlist(selected_pathway[i]),]))),
                  cbind(annot3[annot3$entrezgene %in% unlist(selected_pathway[i]),],sig = rep(-1,nrow(annot3[annot3$entrezgene %in% unlist(selected_pathway[i]),]))))
  index_a_00071 = data.frame(sig = index_a[,4]); rownames(index_a_00071) = index_a$entrezgene
  pathview(gene.data = index_a_00071, pathway.id = pathway_code, species = "bta",kegg.native=T, sign.pos="bottomleft")
  # isthmus
  setwd("/Users/liulihe95/Desktop/metabolomics_0417/")
  dir.create("path_isthmus")
  setwd("/Users/liulihe95/Desktop/metabolomics_0417/path_isthmus")
  index_i = rbind(cbind(annot2[annot2$entrezgene %in% unlist(selected_pathway[i]),],sig = rep(1,nrow(annot2[annot2$entrezgene %in% unlist(selected_pathway[i]),]))),
                  cbind(annot4[annot4$entrezgene %in% unlist(selected_pathway[i]),],sig = rep(-1,nrow(annot4[annot4$entrezgene %in% unlist(selected_pathway[i]),]))))
  index_i_00071 = data.frame(sig = index_i[,4]); rownames(index_i_00071) = index_i$entrezgene
  pathview(gene.data = index_i_00071, pathway.id = c(pathway_code), species = "bta", kegg.native=T, sign.pos="bottomleft")
}
names(selected_pathway)
####################################

i = 2
pathway_code = substr(names(selected_pathway[i]),4,8)
index_i = rbind(cbind(annot2[annot2$entrezgene %in% unlist(selected_pathway[i]),],sig = rep(1,nrow(annot2[annot2$entrezgene %in% unlist(selected_pathway[i]),]))),
                cbind(annot4[annot4$entrezgene %in% unlist(selected_pathway[i]),],sig = rep(-1,nrow(annot4[annot4$entrezgene %in% unlist(selected_pathway[i]),]))))
index_i_00071 = data.frame(sig = index_i[,4]); rownames(index_i_00071) = index_i$entrezgene
pathview.out = pathview(gene.data = index_i_00071, pathway.id = c(pathway_code), species = "bta", kegg.native=T, sign.pos="bottomleft",map.cpdname=TRUE, pdf.size = c(7, 7))
write.csv(index_i_00071,"pathview_webtest.txt")
str(pathview.out)
str(cpd_data_test)
#
cpd_data_test = pathview.out$plot.data.cpd
cpd_data_test$width = 2; cpd_data_test$height = 4
pathview.out2 = pathview(plot.gene.data = index_i_00071, plot.cdp.data =cpd_data_test ,
                        pathway.id = c(pathway_code), species = "bta", kegg.native=T, sign.pos="bottomleft",map.cpdname=TRUE)
str(pathview.out2)
head(pathview.out,6)
getwd()


download.kegg(pathway.id = "00071", species = "bta", kegg.dir = ".",
              file.type=c("xml", "png"))
download.kegg(pathway.id = "01212", species = "bta", kegg.dir = ".",
              file.type=c("xml", "png"))

test_parse =    parseKGML2DataFrame("bta00071.xml",reactions = T)
test_parse_12 = parseKGML2DataFrame("bta01212.xml",reactions = T)

?parseKGML2DataFrame()

str(test_parse)
r??parseKGML2()
?KEGGpathway2Graph2()


index_i_1 = rep(0,length(annot1$entrezgen))
annot2[annot2$entrezgene %in% unlist(selected_pathway[1]),]
annot4[annot4$entrezgene %in% unlist(selected_pathway[1]),]


piel = exp2[exp2$entrezgene %in% sdb$kg.sets$'bta04916 Melanogenesis',]
genes.expresion.diferencial.piel = subset(piel, unos==1)


pathview(gene.data = index_a_input, pathway.id = c("04916","00350"), species = "bta",
         kegg.native=T, sign.pos="bottomleft")

pathview(gene.data = ser, pathway.id = "00350", species = "bta",
         kegg.native=T, sign.pos="bottomleft")

########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################

#View(annot1)
#annot1[table(annot1$ensembl_gene_id) == 2,]
#annot1[annot1$ensembl_gene_id == "ENSBTAG00000027854",]
#annot1 <- annot1[!duplicated(annot1$ensembl_gene_id),]
#setdiff(IDs1,annot1$ensembl_gene_id)
#dim(annot1)


## demo codes
ids1 <- as.vector(row.names(lrt$table))
annot1 <- getBM(attributes=c("ensembl_gene_id", "entrezgene",
                             "external_gene_name"),
                filters="ensembl_gene_id", values=ids1, mart=mart)


length(!is.na(annot1$external_gene_name))

annot1 <- annot1[!duplicated(annot1$ensembl_gene_id),]

exp <- merge(lrt$table, annot1, by.x = 0, by.y = "ensembl_gene_id")
exp1 <- merge(summarized$counts, exp, by.x = 0, by.y = "Row.names")
exp2 <- merge(new.FCC, exp1, by.x = 0, by.y = "Row.names")
exp2 <- exp2[!is.na(exp2$entrezgene),]
exp2 <- exp2[!duplicated(exp2$entrezgene),]
exp2$unos ->FC
exp2$nombres <- rownames(exp2)
as.data.frame(FC)->ser
rownames(ser)<-exp2$entrezgene

#ruta metabolicas
dir.create("path")
setwd("path")
library(pathview)
library(biomaRt)
library(gage)
##obtengo los kegg de bos taurus
kegg.gsets(species = "bta", id.type = "kegg", check.new=FALSE)->sdb
kegg.gs=sdb$kg.sets[sdb$sigmet.id]

#which(IDs1 == overlap_ampulla_LF$ensembl_gene_id)

#IDs1[107];IDs1[199]; IDs1[128]

#which(IDs1 == overlap_ampulla_LF$ensembl_gene_id)

#which(IDs1 == "ENSBTAG00000021905") 


piel = exp2[exp2$entrezgene %in% sdb$kg.sets$'bta04916 Melanogenesis',]
genes.expresion.diferencial.piel = subset(piel, unos==1)


pathview(gene.data = ser, pathway.id = c("04916","00350"), species = "bta",
         kegg.native=T, sign.pos="bottomleft")


pathview(gene.data = ser, pathway.id = "00350", species = "bta",
         kegg.native=T, sign.pos="bottomleft")



#####hago pathway utilizando el FC
exp$logFC ->chu
as.data.frame(chu)->chuu
rownames(chuu)<-rownames(exp)
##melanogenesis
pg <- pathview(gene.data =  chuu, pathway.id = c("04916","00350"),
               species = "bta",
               kegg.native=T, sign.pos="bottomleft", limit =
                 list(gene=max(abs(chuu)), cpd=1))


pathview(gene.data = chuu, pathway.id = "00350", species = "bta",
         kegg.native=T, sign.pos="bottomleft", limit =
           list(gene=max(abs(chuu)), cpd=1))
