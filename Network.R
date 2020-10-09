##### library #####
library(WGCNA);library(ppcor);library(igraph)
library(gdata);library(ggplot2)
library(ggpubr);require(cowplot);library(extrafont);library(dplyr);library(plotly)
library(readxl)
#library(car) ## qqplot but didnt use
library(qqplotr)## qqplot: used 
##################################################################################################
###                                       1. dataprep                                      ######
################################################################################################
setwd("/Users/liulihe95/Desktop/metabolomics_0417")
#install.packages("gdata")
options(stringsAsFactors = FALSE)
networkData = read.xls("metabolomics_data.xlsx",sheet = 1, method=c("tab"), header= T)
rowindex = networkData[,1]
for (i in c(1:(0.5*length(networkData[,1])))){
  rowindex[i] = paste(rowindex[i],i,sep = "_")
  rowindex[i+7] = paste(rowindex[i+7],i,sep = "_")
}
rownames(networkData) = rowindex; networkData = networkData[,-1]
# colnames(networkData)
networkData = data.frame(t(networkData))
datExpr_metab <- as.data.frame(t(networkData));names(datExpr_metab) = rownames(networkData); rownames(datExpr_metab) = rowindex
#str(datExpr_metab)
#dim(datExpr_metab)
#table(colnames(networkData) == rownames(datExpr_metab))
datExpr_low <- datExpr_metab[c(1:7),];datExpr_high <- datExpr_metab[c(8:14),]
dim(datExpr_low);dim(datExpr_high)
###########
### making labels for sub-net --- different chemical properties
allnames_label = colnames(datExpr_low)
label1 = c(1:26) # Amino acids and Biogenic Amines 
label2= c(27:66) # Acylcarnitines 
label3 = c(67:80)  # Lysophosphatidylcholines
label4 = c(81:153)  # Phosphatidylcholines
label5 = c(154:167)  # Sphingomyelins 
label6 = c(168)    # Hexoses 
label7 = c(167:179) # Eicosanoids & oxidation products of polyunsaturated fatty acids
length(label7)

##################################################################
#### Function pre ###############################################
##################################################################
check_zero = function(datExpr1,datExpr2,start,stop,thres){
  label_low1 = as.numeric() ; label_low2 = as.character()
  label_high1 = as.numeric()  ; label_high2 = as.character()
  for (i in c(1:length(allnames_label))){
    test1 = length(which(datExpr1[,i]==0))
    test2 = length(which(datExpr2[,i]==0))
    if (test1  >= 7) {
      label_low1[i] = i
      label_low2[i] = allnames_label[i]
    }
    else if (test2 >= 7){
      label_high1[i] = i
      label_high2[i] = allnames_label[i]
    }
    label_low1 = label_low1[!is.na(label_low1)]
    label_low2 = label_low2[!is.na(label_low2)]
    label_high1 = label_high1[!is.na(label_high1)]
    label_high2 = label_high2[!is.na(label_high2)]
    remove = unique(c(label_low1,label_high1))
  }
  Results = list(position1 = label_low1, position2 = label_low2, name1 = label_high1, name2 = label_high2,may_remove = remove)
  return(Results)
}
checkzero_result = check_zero(datExpr_low,datExpr_high,1,length(datExpr_high),7)
remove_index = checkzero_result$may_remove; remove_index## remove 9 and 171, now we have 179 - 2 = 177 nodes
rownames(networkData)[remove_index]
##################################################################
get_NC_cor = function(datExpr1,datExpr2,r_thres,p_thres){
  #calcu of matrx1
  cormatr1 <- cor(datExpr1)
  adjmatr1 = matrix(1,ncol(datExpr1),ncol(datExpr1))
  colnames(adjmatr1) = colnames(datExpr1)
  rownames(adjmatr1) = colnames(datExpr1)
  adjmatr1[abs(cormatr1) < r_thres] = 0
  #calcu of matrx2
  cormatr2 <- cor(datExpr2)
  adjmatr2 = matrix(1,ncol(datExpr2),ncol(datExpr2))
  colnames(adjmatr2) = colnames(datExpr2)
  rownames(adjmatr2) = colnames(datExpr2)
  adjmatr2[abs(cormatr2) < r_thres] = 0
  #use threshold and get pvalue matrix
  pvalmatr1 = matrix(0,ncol(datExpr1),ncol(datExpr1))
  pvalmatr2 = matrix(0,ncol(datExpr2),ncol(datExpr2))
  for(i in 1:(ncol(datExpr1)-1)){
    for (j in c((i+1):ncol(datExpr1))){
      r1 = cor.test(datExpr1[,i],datExpr1[,j])
      r2 = cor.test(datExpr2[,i],datExpr2[,j])
      pvalmatr1[i,j] = pvalmatr1[j,i] = r1$p.value
      pvalmatr2[i,j] = pvalmatr2[j,i] = r2$p.value
      if(r1$p.value >= p_thres){adjmatr1[i,j] = adjmatr1[j,i] = 0}
      if(r2$p.value >= p_thres){adjmatr2[i,j] = adjmatr2[j,i] = 0}
    }
  }
  #get all the basic NC
  NC1 = conformityBasedNetworkConcepts(adjmatr1)
  NC2 = conformityBasedNetworkConcepts(adjmatr2)
  #combine, rank and show
  basic_results = data.frame(density1 = NC1$fundamentalNCs$Density,
                             density2 = NC2$fundamentalNCs$Density,
                             centralization1 = NC1$fundamentalNCs$Centralization,
                             centralization2 = NC2$fundamentalNCs$Centralization,
                             heterogeneity1 = NC1$fundamentalNCs$Heterogeneity,
                             heterogeneity2 = NC2$fundamentalNCs$Heterogeneity)
  change_results = data.frame(#gene = colnames(datExpr1),
    con_1 = NC1$fundamentalNCs$Connectivity,
    scl_con_1 = NC1$fundamentalNCs$ScaledConnectivity,
    con_2 = NC2$fundamentalNCs$Connectivity,
    scl_con_2 = NC2$fundamentalNCs$ScaledConnectivity,
    con_change = -(NC1$fundamentalNCs$Connectivity - NC2$fundamentalNCs$Connectivity),
    scl_con_change = -(NC1$fundamentalNCs$ScaledConnectivity - NC2$fundamentalNCs$ScaledConnectivity),
    rank_scl_con = rank(-abs(NC1$fundamentalNCs$ScaledConnectivity-NC2$fundamentalNCs$ScaledConnectivity)),
    cls_coef_1 = NC1$fundamentalNCs$ClusterCoef,
    cls_coef_2 = NC2$fundamentalNCs$ClusterCoef,
    clst_coef_change = c(NC1$fundamentalNCs$ClusterCoef - NC2$fundamentalNCs$ClusterCoef),
    rank_clstcoef = rank(-abs(NC1$fundamentalNCs$ClusterCoef-NC2$fundamentalNCs$ClusterCoef)))
  Results = list(NC1 = NC1, NC2 = NC2, cormatr_ref = cormatr1, cormatr_test = cormatr2, pvalue1 = pvalmatr1,pvalue2 = pvalmatr2 ,adjmatr_ref = adjmatr1, adjmatr_test = adjmatr2, basic = basic_results,change = change_results)
  return(Results)
}


##################################################################
######                    check for missing value                  ##########
# gsg_all_sperm = goodSamplesGenes(datExpr_metab, verbose = 3);
# gsg_all_sperm$allOK
####   distance between samples / outliers
# sample_Tree_low = hclust(dist(datExpr_low), method = "average")
# sample_Tree_high = hclust(dist(datExpr_high), method = "average")
# plot(sample_Tree_low, main = "Sample clustering to detect outliers low", sub="", xlab="", cex.lab = 1.5, 
#      cex.axis = 1.5, cex.main = 2)
# plot(sample_Tree_high, main = "Sample clustering to detect outliers high", sub="", xlab="", cex.lab = 1.5, 
#      cex.axis = 1.5, cex.main = 2)
##################################################################
datExpr_low_new = datExpr_low[,-remove_index]
datExpr_high_new = datExpr_high[,-remove_index]

################### Get network #################################
Results_metab_cor = get_NC_cor(datExpr_low_new,datExpr_high_new,0.5,0.05)
#label3 #Results_metab_cor_testsub3 = get_NC_cor(datExpr_low_new[,label3],datExpr_high_new[,label3],0.5,0.05)
# spearman, (with ties (no exact p value))
#Results_metab_spearman = get_NC_spearman(datExpr_low_new,datExpr_high_new,0.5,0.05)
# Results_metab_cor
Results_metab_cor$basic
length(Results_metab_cor$NC1$fundamentalNCs$Connectivity)

###########################################################################
########                   2. betwork statistic                  ############
##########################################################################
### 2.connectivity and cluster coef 
######  connectivity (of each node)  vector  ##################
Con_ref= rowSums(Results_metab_cor$adjmatr_ref) - 1 
Con_test= rowSums(Results_metab_cor$adjmatr_test) - 1
hist(Con_ref)
hist(Con_test)
dev.off()
###### mean connectivity  ##############
meanCon_ref = sum(Con_ref)/ncol(datExpr_low_new);meanCon_ref
meanCon_test = sum(Con_test)/ncol(datExpr_high_new);meanCon_test

######     density      #################
density_ref = sum(vectorizeMatrix(Results_metab_cor$adjmatr_ref))/(0.5*ncol(datExpr_low_new)*(ncol(datExpr_low)-1));density_ref
density_test = sum(vectorizeMatrix(Results_metab_cor$adjmatr_test))/(0.5*ncol(datExpr_high_new)*(ncol(datExpr_high_new)-1));density_test
######  coef (of each node)  vector  ##################
clstcoef_ref= Results_metab_cor$NC1$fundamentalNCs$ClusterCoef
clstcoef_test= Results_metab_cor$NC2$fundamentalNCs$ClusterCoef
###### mean coef  ##################
meanClstcoef_ref = sum(clstcoef_ref)/ncol(datExpr_low_new);meanClstcoef_ref
meanClstcoef_test = sum(clstcoef_test)/ncol(datExpr_high_new);meanClstcoef_test

### 3.top 10 of connectivity/clustercoeffcient ###
########## assemble dataset 1 _ con
# topgene_metab_con = data.frame(
#   Con_ref = Results_metab_cor$change$con_1,
#   Con_rank_ref = rank(-Results_metab_cor$change$con_1,ties.method = "min"),
#   #Con_ref_scl = Results_sperm_cor$change$scl_con_1,
#   Con_test = Results_metab_cor$change$con_2,
#   #Con_ref_scl = Results_sperm_cor$change$scl_con_1,
#   Con_rank_test = rank(-Results_metab_cor$change$con_2,ties.method = "min"),
#   Con_Change = Results_metab_cor$change$con_change,
#   ConChange_rank =rank(-(abs(Results_metab_cor$change$con_1 - Results_metab_cor$change$con_2)),ties.method = "min")
# )
# rownames(topgene_metab_con) = colnames(datExpr_high_new)
# ######### assemble dataset 2 _ clscoef
# topgene_metab_clscoef = data.frame(
#   Clscoef_ref = Results_metab_cor$change$cls_coef_1 ,
#   Clscoef_rank_ref = rank(-Results_metab_cor$change$cls_coef_1,ties.method = "min"),
#   Clscoef_test = Results_metab_cor$change$cls_coef_2,
#   Clscoef_rank_test = rank(-Results_metab_cor$change$cls_coef_2,ties.method = "min"),
#   Clscoef_Change = Results_metab_cor$change$clst_coef_change,
#   ConChange_rank =  rank(-(abs(Results_metab_cor$change$cls_coef_1 - Results_metab_cor$change$cls_coef_2)),ties.method = "min")
# )
# rownames(topgene_metab_clscoef) = colnames(datExpr_high_new)
# ###################################################################################################
# ## proportion analysis - how many common neighbors  ### ???
# Results_metab_cor$change[,1:2]
# sig_prop = which(rank(Results_metab_cor$change$proportion,ties.method = "min") <= 65)
# sig_conchange = which(rank(-abs(Results_metab_cor$change$con_change),ties.method = "min") <= 65)
# sig_clstchange = which(rank(-abs(Results_metab_cor$change$clst_coef_change),ties.method = "min") <= 65)
# summary(Results_metab_cor$change$proportion)
# #n
# length(which(Results_metab_cor$change$proportion == 0))
# 
# all_three_low = intersect(intersect(sig_prop, sig_conchange),sig_clstchange)
# rownames(test_ranking)[all_three_low]
# 
# ###################################################################################################
# ## rewiring score analysis
# adj1 = Results_metab_cor$adjmatr_ref; adj2 = Results_metab_cor$adjmatr_test
# # mean(nonzero(1,2,0))
# 
# ##############################################################################################################################
# # Ranking test
# prop_all = Results_metab_cor$change$proportion
# conchange_all = abs(Results_metab_cor$change$con_change)
# clstchange_all = abs(Results_metab_cor$change$clst_coef_change)
# 
# test_ranking <- do.call('cbind', list(prop_all,conchange_all,clstchange_all))
# rownames(test_ranking) = rownames(Results_metab_cor$change)
# ranktest = rankPvalue(test_ranking,ties.method = "min")
# names(ranktest)
# ranktest[ranktest$qValueLowRank<=0.1,]
# ranktest[ranktest$pValueLowRank<=0.05,]
# 
# 
# cor(Results_metab_cor$change$con_change,Results_metab_cor$change$clst_coef_change,method = "spearman")
# cor(Results_metab_cor$change$con_change,Results_metab_cor$change$proportion,method = "spearman")
# cor(Results_metab_cor$change$clst_coef_change,Results_metab_cor$change$proportion,method = "spearman")
# 


# values in the columns of datS are independent. 
# This allows us to find objects (rows) with consistently high (or low) values across the columns.
# the function rankPvalue can be used to calculate a p-value for this occurrence.
# central limit theorem (referred to as percentile rank method)

# defined as the sum of the percentile ranks

# 
# #### define function 
# SelectGene_un_cor = function(dataset,topnumber){
#   index1 = dataset[,2] %in% c(1:topnumber)
#   index2 = dataset[,4] %in% c(1:topnumber)
#   index3 = dataset[,6] %in% c(1:topnumber)
#   summary1 = dataset[index1,];summary1 = summary1[order(summary1[,2]),]
#   summary2 = dataset[index2,];summary2 = summary2[order(summary2[,4]),]
#   summary3 = dataset[index3,];summary3 = summary3[order(summary3[,6]),]
#   summary = list(
#     ref = summary1,
#     test = summary2,
#     change =summary3
#   )
#   return(summary)
# }
# ##### example result, both in top 22 - con and clustercoef
# top10_con_change = rownames(SelectGene_un_cor(topgene_metab_con,22)[[3]])
# top10_ClusterCoef_change = rownames(SelectGene_un_cor(topgene_metab_clscoef,22)[[3]])
# intersect(top10_con_change,top10_ClusterCoef_change)
# Results_metab_cor$basic

###################################################################################################
####                                       3. Plotting                                      ######
#################################################################################################
#################################        1. generating dataset ####################
#####  con --- test_combine_dataset
ref = data.frame(
  connectivity = as.numeric(Results_metab_cor$NC1$fundamentalNCs$Connectivity),
  category = rep("LF-LCL",ncol(datExpr_low_new)))
test = data.frame(
  connectivity = as.numeric(Results_metab_cor$NC2$fundamentalNCs$Connectivity),
  category = rep("SF-SCL",ncol(datExpr_high_new)))
test_combine_dataset <- do.call('rbind', list(ref,test))
str(test_combine_dataset)
# table(test_combine_dataset$category)

### clst coef --- test_combine_dataset_clstcoef
ref_clscoef = data.frame(
  clstcoef = as.numeric(Results_metab_cor$NC1$fundamentalNCs$ClusterCoef),
  category = rep("LF-LCL",ncol(datExpr_low_new)))
test_clscoef = data.frame(
  clstcoef = as.numeric(Results_metab_cor$NC2$fundamentalNCs$ClusterCoef),
  category = rep("SF-SCL",ncol(datExpr_high_new)))
test_combine_dataset_clstcoef <- do.call('rbind', list(ref_clscoef,test_clscoef))
str(test_combine_dataset_clstcoef)

Results_metab_cor$basic

library(extrafont)
##################################        2. plotting         ##########################
#plot3
#max(test_combine_dataset$connectivity)
plot3 = ggplot(test_combine_dataset, aes(x=connectivity, fill=category)) +
  geom_histogram(binwidth=1,alpha=0.75, position="identity", aes(y = ..count..), color="black") +
  # geom_density(alpha=0.6,trim = F) +
  xlim(-2,(max(test_combine_dataset$connectivity)+5))+
  geom_vline(aes(xintercept=meanCon_ref), color="black", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=meanCon_test), color="blue", linetype="dashed", size=1) +
  theme_gray()+
  #scale_fill_discrete(name="Groups")+
  theme(legend.position="top",legend.title = element_text(size = 10, color = "black",face = "bold",vjust = 0.5, hjust = 0.5))+ #family = "Microsoft Sans Serif",
  scale_fill_manual(name="Groups",values = c("blue", "red")) +
  #theme(legend.position="top")+
  #scale_fill_manual(name="Experimental\nCondition",breaks=c("ref", "test"),labels=c("Small","Large"))+
  ##scale_fill_manual(name="Groups",values=c("blue", "red"))+#,# name="Experimental\nCondition",breaks=c("ctrl", "trt1", "trt2"),labels=c("Control", "Treatment 1", "Treatment 2"))
  labs(title="Distribution of Connectivity", x="Connectivity", y = "Frequency")+
  theme(axis.text.x = element_text(size = 15, family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 15,family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x = element_text(size = 20,family = "Microsoft Sans Serif",color = "black",vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 20, color = "black",family = "Microsoft Sans Serif", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(size = 20, family = "Microsoft Sans Serif",color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))
plot3
#  plot4 
plot4 = ggplot(test_combine_dataset_clstcoef, aes(x=clstcoef, fill=category)) +
  geom_histogram(binwidth=.01,alpha=0.75, position="identity", aes(y = ..count..), color="black") +
  #  geom_density(alpha=0.6,trim = F) +
  xlim(-0.05,1.05)+
  geom_vline(aes(xintercept=meanClstcoef_ref), color="black", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=meanClstcoef_test), color="blue", linetype="dashed", size=1) +
  theme_gray()+
  #scale_fill_discrete(name="Groups")+
  scale_fill_manual(name="Groups",values = c("blue", "red")) +
  #theme(legend.position="None")+
  theme(legend.position="top",legend.title = element_text(size = 10, color = "black", face = "bold",vjust = 0.5, hjust = 0.5))+ #,family = "Microsoft Sans Serif"
  labs(title="Distribution of Cluster Coefficient", x="Cluster Coefficient", y = "Frequency")+
  theme(axis.text.x = element_text(size = 15, family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 15,family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x = element_text(size = 20,family = "Microsoft Sans Serif",color = "black",vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 20, color = "black",family = "Microsoft Sans Serif", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(size = 20, family = "Microsoft Sans Serif",color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))
plot4

dev.off()

tiff("Metab_network_plot.tiff", width = 14, height = 12, units = 'in', res = 300)
plot_grid(plot3, plot4, align = c("v"),labels = c("A","B"), nrow = 2,label_size= 20, label_colour = "darkgreen")
#tiff("presentation_figure.tiff", width = 14, height = 12, units = 'in', res = 300)
#plot_grid(plot3,plot4, align = c("v"),labels = c("A","B"), nrow = 2,label_size= 20, label_colour = "darkgreen")
dev.off()

###################################################################################################
####                                       3. pre-difine (chemical category)                ######
#################################################################################################
label1 = c(1:26);label1 = label1[-intersect(label1,remove_index)] # Amino acids and Biogenic Amines 
label2 = c(27:66);#label2 = label2[-intersect(label2,remove_index)]; # Acylcarnitines 
label3 = c(67:80);#label3 = label3[-intersect(label3,remove_index)];  # Lysophosphatidylcholines
label4 = c(81:153);#label4 = label4[-intersect(label4,remove_index)];  # Phosphatidylcholines
label5 = c(154:167);#label5 = label5[-intersect(label5,remove_index)]; # Sphingomyelins 
# label6 = c(168)    # Hexoses 
label7 = c(167:179); label7 = label7[-intersect(label7,remove_index)]; # Amino acids and Biogenic Amines 
label7 = c(167:170,172:179); #label7 = label7[-171] # Amino acids and Biogenic Amines 
# 
sub1_names = colnames(datExpr_high)[label1]
sub2_names = colnames(datExpr_high)[label2]
sub3_names = colnames(datExpr_high)[label3]
sub4_names = colnames(datExpr_high)[label4]
sub5_names = colnames(datExpr_high)[label5]
sub7_names = colnames(datExpr_high)[label7]
# 
#Results_metab_cor_sub1$basic
#Results_metab_cor_sub1$change
#
#Results_metab_cor$change[sub3_names,]
# remove_index
### sub - datasets
###################################################################################################
####                              4. subnet of pre-defined module                           ######
#################################################################################################
########## 
allnames_label[label1]
label1 = c(1:26);label1 = label1[-intersect(label1,remove_index)] # Amino acids and Biogenic Amines 
label2 = c(27:66);#label2 = label2[-intersect(label2,remove_index)]; # Acylcarnitines 
label3 = c(67:80);#label3 = label3[-intersect(label3,remove_index)];  # Lysophosphatidylcholines
label4 = c(81:153);#label4 = label4[-intersect(label4,remove_index)];  # Phosphatidylcholines
label5 = c(154:167);#label5 = label5[-intersect(label5,remove_index)]; # Sphingomyelins 
# label6 = c(168)    # Hexoses 
label7 = c(167:170,172:179); #label7 = label7[-171] # Amino acids and Biogenic Amines 
## dataset pre
test_data_sub1_low = datExpr_low[,label1]; test_data_sub1_high = datExpr_high[,label1]
test_data_sub2_low = datExpr_low[,label2]; test_data_sub2_high = datExpr_high[,label2]
test_data_sub3_low = datExpr_low[,label3]; test_data_sub3_high = datExpr_high[,label3]
test_data_sub4_low = datExpr_low[,label4]; test_data_sub4_high = datExpr_high[,label4]
test_data_sub5_low = datExpr_low[,label5]; test_data_sub5_high = datExpr_high[,label5]
test_data_sub7_low = datExpr_low[,label7]; test_data_sub7_high = datExpr_high[,label7]

### results pre
testResults_cor_sub1 = get_NC_cor(test_data_sub1_low,test_data_sub1_high,0.5,0.05)
testResults_spearmen_sub1 = get_NC_spearman(test_data_sub1_low,test_data_sub1_high,0.5,0.05)
testResults_cor_sub2 = get_NC_cor(test_data_sub2_low,test_data_sub2_high,0.5,0.05)
testResults_spearmen_sub2 = get_NC_spearman(test_data_sub2_low,test_data_sub2_high,0.5,0.05)
#
testResults_cor_sub3 = get_NC_cor_modf(test_data_sub3_low,test_data_sub3_high,0.5,0.05)
testResults_spearmen_sub3 = get_NC_spearman(test_data_sub3_low,test_data_sub3_high,0.5,0.05)
testResults_cor_sub4 = get_NC_cor(test_data_sub4_low,test_data_sub4_high,0.5,0.05)
testResults_spearmen_sub4 = get_NC_spearman(test_data_sub4_low,test_data_sub4_high,0.5,0.05)
#
testResults_cor_sub5 = get_NC_cor(test_data_sub5_low,test_data_sub5_high,0.5,0.05)
testResults_spearmen_sub5 = get_NC_spearman(test_data_sub5_low,test_data_sub5_high,0.5,0.05)
testResults_cor_sub7 = get_NC_cor(test_data_sub7_low,test_data_sub7_high,0.5,0.05)
testResults_spearmen_sub7 = get_NC_spearman_modf(datExpr_low[,label7],test_data_sub7_high,0.5,0.05)



# plotting - all 
plot3_1 = ggplot(test_combine_dataset, aes(x=connectivity, fill=category)) +
  geom_histogram(binwidth=1,alpha=0.6, position="identity", aes(y = ..count..), color="black") +
  # geom_density(alpha=0.6,trim = F) +
  #xlim(0,(max(test_combine_dataset$connectivity)+5))+
  geom_vline(aes(xintercept=meanCon_ref), color="black", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=meanCon_test), color="blue", linetype="dashed", size=1) +
  theme_gray()+
  theme(legend.position="None")+
  labs(title="all_metab_con_pearson", x="Connectivity", y = "Frequency")+
  theme(axis.text.x = element_text(size = 15, family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 15,family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x = element_text(size = 15,family = "Microsoft Sans Serif",color = "black",vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 15, color = "black",family = "Microsoft Sans Serif", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(size = 20, family = "Microsoft Sans Serif",color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))
plot4_1  = ggplot(test_combine_dataset_clstcoef, aes(x=clstcoef, fill=category)) +
  geom_histogram(binwidth=.01,alpha=0.6, position="identity", aes(y = ..count..), color="black") +
  #  geom_density(alpha=0.6,trim = F) +
  xlim(-0.1,1.1)+
  geom_vline(aes(xintercept=meanClstcoef_ref), color="black", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=meanClstcoef_test), color="blue", linetype="dashed", size=1) +
  theme_gray()+
  theme(legend.position="None")+
  labs(title="all_metab_clst_pearson", x="Cluster Coefficient", y = "Frequency")+
  theme(axis.text.x = element_text(size = 15, family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 15,family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x = element_text(size = 15,family = "Microsoft Sans Serif",color = "black",vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 15, color = "black",family = "Microsoft Sans Serif", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(size = 20, family = "Microsoft Sans Serif",color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))
plot3_spearman = Plot_compare_distrib_con(Results_metab_spearman$change$con_1,Results_metab_spearman$change$con_2,"all_metab_con_spearman")
plot4_spearman = Plot_compare_distrib_clst(Results_metab_spearman$change$cls_coef_1,Results_metab_spearman$change$cls_coef_2,"all_metab_clst_spearman")

multiplot(plot3_1,plot4_1,plot3_spearman,plot4_spearman,cols = 2)
dev.off()
#  plots pre
label1 = c(1:26) # Amino acids and Biogenic Amines 
label2 = c(27:66) # Acylcarnitines 
label3 = c(67:80)  # Lysophosphatidylcholines
label4 = c(81:153)# Phosphatidylcholines
label5 = c(154:167)  # Sphingomyelins 
label6 = c(168)    # Hexoses 
label7 = c(167:179) # Eicosanoids & oxidation products of polyunsaturated fatty acids

######################################################### cor #########################################################
test_cor_con1 = Plot_compare_distrib_con(testResults_cor_sub1$change$con_1,testResults_cor_sub1$change$con_2,"sub1_Amino acids and Biogenic Amines_con")
test_cor_clst1 = Plot_compare_distrib_clst(testResults_cor_sub1$change$cls_coef_1,testResults_cor_sub1$change$cls_coef_2,"sub1_Amino acids and Biogenic Amines_clst")
#sub2
test_cor_con2 = Plot_compare_distrib_con(testResults_cor_sub2$change$con_1,testResults_cor_sub2$change$con_2,"sub2_Acylcarnitines_con")
test_cor_clst2 = Plot_compare_distrib_clst(testResults_cor_sub2$change$cls_coef_1,testResults_cor_sub2$change$cls_coef_2,"sub2_Acylcarnitines_clst")
tiff("subplot1.tiff", width = 14, height = 12, units = 'in', res = 300)
multiplot(test_cor_con1,test_cor_clst1,test_cor_con2,test_cor_clst2,cols = 2)
dev.off()
# Plot_compare_distrib_clst(Results_metab_cor$change$cls_coef_1,Results_metab_cor$change$cls_coef_2,"sub2_clst")
#sub3
test_cor_con3 = Plot_compare_distrib_con(testResults_cor_sub3$change$con_1,testResults_cor_sub3$change$con_2,"sub3_Lysophosphatidylcholines_con")
test_cor_clst3 = Plot_compare_distrib_clst(testResults_cor_sub3$change$cls_coef_1,testResults_cor_sub3$change$cls_coef_2,"sub3_Lysophosphatidylcholines_clst")
#sub4
test_cor_con4 = Plot_compare_distrib_con(testResults_cor_sub4$change$con_1,testResults_cor_sub4$change$con_2,"sub4_Phosphatidylcholines_con")
test_cor_clst4 = Plot_compare_distrib_clst(testResults_cor_sub4$change$cls_coef_1,testResults_cor_sub4$change$cls_coef_2,"sub4_Phosphatidylcholines_clst")
tiff("subplot2.tiff", width = 14, height = 12, units = 'in', res = 300)
multiplot(test_cor_con3,test_cor_clst3,test_cor_con4,test_cor_clst4,cols = 2)
dev.off()
#sub5
test_cor_con5 = Plot_compare_distrib_con(testResults_cor_sub5$change$con_1,testResults_cor_sub5$change$con_2,"sub5_Sphingomyelins_con")
test_cor_clst5 = Plot_compare_distrib_clst(testResults_cor_sub5$change$cls_coef_1,testResults_cor_sub5$change$cls_coef_2,"sub5_Sphingomyelins_clst")
#sub7
test_cor_con7 = Plot_compare_distrib_con(testResults_cor_sub7$change$con_1,testResults_cor_sub7$change$con_2,"sub7_Eicosanoids & ps-fa_con")
test_cor_clst7 = Plot_compare_distrib_clst(testResults_cor_sub7$change$cls_coef_1,testResults_cor_sub7$change$cls_coef_2,"sub7_Eicosanoids & ps-fa_clst")
tiff("subplot3.tiff", width = 14, height = 12, units = 'in', res = 300)
multiplot(test_cor_con5,test_cor_clst5,test_cor_con7,test_cor_clst7,cols = 2)
dev.off()


######################################################### spearman  #########################################################
# sub1
test_spearmen_con1 = Plot_compare_distrib_con(testResults_spearmen_sub1$change$con_1,testResults_spearmen_sub1$change$con_2,"sub1_Amino acids and Biogenic Amines_con")
test_spearmen_clst1 = Plot_compare_distrib_clst(testResults_spearmen_sub1$change$cls_coef_1,testResults_spearmen_sub1$change$cls_coef_2,"sub1_Amino acids and Biogenic Amines_clst")
#sub2
test_spearmen_con2 = Plot_compare_distrib_con(testResults_spearmen_sub2$change$con_1,testResults_spearmen_sub2$change$con_2,"sub2_Acylcarnitines_con")
test_spearmen_clst2 = Plot_compare_distrib_clst(testResults_spearmen_sub2$change$cls_coef_1,testResults_spearmen_sub2$change$cls_coef_2,"sub2_Acylcarnitines_clst")
multiplot(test_spearmen_con1,test_spearmen_clst1,test_spearmen_con2,test_spearmen_clst2,cols = 2)
#sub3
test_spearmen_con3 = Plot_compare_distrib_con(testResults_spearmen_sub3$change$con_1,testResults_spearmen_sub3$change$con_2,"sub3_Lysophosphatidylcholines_con")
test_spearmen_clst3 = Plot_compare_distrib_clst(testResults_spearmen_sub3$change$cls_coef_1,testResults_spearmen_sub3$change$cls_coef_2,"sub3_Lysophosphatidylcholines_clst")
#sub4
test_spearmen_con4 = Plot_compare_distrib_con(testResults_spearmen_sub4$change$con_1,testResults_spearmen_sub4$change$con_2,"sub4_Phosphatidylcholines_con")
test_spearmen_clst4 = Plot_compare_distrib_clst(testResults_spearmen_sub4$change$cls_coef_1,testResults_spearmen_sub4$change$cls_coef_2,"sub4_Phosphatidylcholines_clst")
multiplot(test_spearmen_con3,test_spearmen_clst3,test_spearmen_con4,test_spearmen_clst4,cols = 2)
#sub5
test_spearmen_con5 = Plot_compare_distrib_con(testResults_spearmen_sub5$change$con_1,testResults_spearmen_sub5$change$con_2,"sub5_Sphingomyelins_con")
test_spearmen_clst5 = Plot_compare_distrib_clst(testResults_spearmen_sub5$change$cls_coef_1,testResults_spearmen_sub5$change$cls_coef_2,"sub5_Sphingomyelins_clst")
mean(testResults_spearmen_sub5$change$cls_coef_1)
#sub7
test_spearmen_con7 = Plot_compare_distrib_con(testResults_spearmen_sub7$change$con_1,testResults_spearmen_sub7$change$con_2,"sub7_Eicosanoids & ps-fa_con")
test_spearmen_clst7 = Plot_compare_distrib_clst(testResults_spearmen_sub7$change$cls_coef_1,testResults_spearmen_sub7$change$cls_coef_2,"sub7_Eicosanoids & ps-fa_clst")
multiplot(test_spearmen_con5,test_spearmen_clst5,test_spearmen_con7,test_spearmen_clst7,cols = 2)
dev.off()
####

###################################################################################################
####                                       5. Expo to cyto                                  ######
#################################################################################################
adjmatr_low = Results_metab_cor$adjmatr_ref
adjmatr_high = Results_metab_cor$adjmatr_test
diag(adjmatr_low) = 0;diag(adjmatr_high) = 0
#
Net_low = graph.adjacency(adjmatr_low,mode = "undirected", weighted = NULL)
Net_high = graph.adjacency(adjmatr_high,mode = "undirected", weighted = NULL)
#
write_graph(Net_low, "Net_ref_metab.text", format = "edgelist")
write_graph(Net_high, "Net_test_metab.text", format = "edgelist")
#
Net_low_list = read.csv("Net_ref_metab.text",sep = "",col.names = c("source","end"))
Net_high_list = read.csv("Net_test_metab.text",sep = "",col.names = c("source","end"))

Net_low_list_m = Net_low_list + 1
Net_high_list_m = Net_high_list +1

# test the if the number and the names match
# length(which(Net_low_list_m[,1] == "85")) + length(which(Net_low_list_m[,2] == "85"))

# get all names (remove 2)
allnames_label_new = allnames_label[-remove_index]
# do massage of dataset / match the node
# get order accroding to "con1" in the result_compilation
names_with_order = rownames(Results_metab_cor$change[,1:2][order(-Results_metab_cor$change[,1]),])



order_assign_x = cbind(allnames_label_new.x = names_with_order,factor(seq(length(names_with_order),1,by= -1),ordered = F))
order_assign_y = cbind(allnames_label_new.y = names_with_order,factor(seq(length(names_with_order),1,by= -1),ordered = F))

#identical(order_assign_x[,2],order_assign_y[,2])

#order_assign_x[which(order_assign_x[,1] =="Leu"),]
#step5_2[step5_2$allnames_label_new.y=="Leu",]
#step5_2[step5_2$allnames_label_new.x=="Leu",]

# match and assign names with the Gene ID
allID_sub_annot_1 = cbind(allnames_label_new, source = factor(seq(1,length(allnames_label_new),by= 1),ordered = F)) # for the use of name "source"
allID_sub_annot_2 = cbind(allnames_label_new, end = factor(seq(1,length(allnames_label_new),by= 1),ordered = F)) # for the use of name "end"
allID_sub_annot = cbind(allnames_label_new, index = factor(seq(1,length(allnames_label_new),by= 1),ordered = F))


Match_Name2ID = data.frame(allID_sub_annot_1)
names(Match_Name2ID) = c('Names','ID')
Match_Name2ID$ID = as.integer(Match_Name2ID$ID)

Assign_plot_order = data.frame(order_assign_x)
names(Assign_plot_order) = c('Names','Order')


library(tidyverse)
edgelist_large = Net_low_list_m %>% as_tibble() %>% 
  dplyr::left_join(Match_Name2ID, by = c('source'='ID')) %>% 
  dplyr::rename(Names_x = Names) %>%
  dplyr::left_join(Match_Name2ID, by = c('end'='ID')) %>% 
  dplyr::rename(Names_y = Names) %>% 
  dplyr::left_join(Assign_plot_order, by = c('Names_x'='Names'))

edgelist_small = Net_high_list_m %>% as_tibble() %>% 
  dplyr::left_join(Match_Name2ID, by = c('source'='ID')) %>% 
  dplyr::rename(Names_x = Names) %>%
  dplyr::left_join(Match_Name2ID, by = c('end'='ID')) %>% 
  dplyr::rename(Names_y = Names) %>% 
  dplyr::left_join(Assign_plot_order, by = c('Names_x'='Names'))

write.csv(edgelist_large,"Final_edgelist_cyto_large.txt",quote = F)
write.csv(edgelist_small,"Final_edgelist_cyto_small.txt",quote = F)


Assign_plot_order %>% 
  dplyr::filter(Names == 'AA')




# # low
# step1 = merge(allID_sub_annot_1, Net_low_list_m, by ="source")
# step2 = merge(allID_sub_annot_2, step1, by ="end");
# # high
# step3 = merge(allID_sub_annot_1, Net_high_list_m, by ="source")
# step4 = merge(allID_sub_annot_2, step3, by ="end")
# 
# # check
# # length(which(Net_low_list[,1] == "85")) + length(which(Net_low_list[,2] == "85"))
# # length(which(step2[,1] == "85")) + length(which(step2[,3] == "85"))
# step5 = merge(order_assign_x,step2,by ="allnames_label_new.x");
# step5_2 = merge(order_assign_y, step5, by ="allnames_label_new.y");
# 
# step5_2[step5_2$allnames_label_new.y=="PC.aa.C32.2",]
# step5_2[step5_2$allnames_label_new.x=="PC.aa.C32.2",]
# 
# 
# colnames(step5_2)[colnames(step5_2)=="V2.x"] <- "plot_order_x";colnames(step5_2)[colnames(step5_2)=="V2.y"] <- "plot_order_y"
# step5_2[,"plot_order_x"] = as.numeric(step5_2[,"plot_order_x"] );step5_2[,"plot_order_y"] = as.numeric(step5_2[,"plot_order_y"] )
# str(step5_2)
# length(unique(step5_2$allnames_label_new.y))
# length(unique(step5_2$allnames_label_new.x))
# 
# # check
# # length(which(step5_2[,1] == "PC.aa.C32.0")) + length(which(step5_2[,3] == "PC.aa.C32.0"))
# 
# step6 = merge(order_assign_x, step4, by ="allnames_label_new.x");
# step6_2 = merge(order_assign_y, step6, by ="allnames_label_new.y");
# colnames(step6_2)[colnames(step6_2)=="V2.x"] <- "plot_order_x";colnames(step6_2)[colnames(step6_2)=="V2.y"] <- "plot_order_y"
# step6_2[,"plot_order_x"] = as.numeric(step6_2[,"plot_order_x"] );step6_2[,"plot_order_y"] = as.numeric(step6_2[,"plot_order_y"] )
# str(step6_2)
# # check
# # length(which(step6_2[,1] == "PC.aa.C32.0")) + length(which(step6_2[,3] == "PC.aa.C32.0"))
# Final_edgelist_low = step5_2;Final_edgelist_high = step6_2
# 
# colnames(Final_edgelist_low)
# Final_edgelist_low_cyto = Final_edgelist_low[,c(1:4)]
# head(Final_edgelist_low_cyto,2)
# 
# # str(Final_edgelist_low)
# Final_edgelist_low_cyto[Final_edgelist_low$allnames_label_new.x == "Leu",]
# Final_edgelist_low_cyto[Final_edgelist_low$allnames_label_new.y == "Leu",]
# 
# 
# # check
# # dim(step5_2);dim(Net_low_list)
# # length(which(Final_edgelist_low[,1] == "PC.aa.C32.0")) + length(which(Final_edgelist_low[,3] == "PC.aa.C32.0"))
# # length(which(Net_low_list[,1] == "85")) + length(which(Net_low_list[,2] == "85"))
# # str(step5_2)
# # str(Net_low_list)
# 
# write.csv(Final_edgelist_low,"Final_edgelist_cyto_low.txt",quote = F)
# write.csv(Final_edgelist_low_cyto,"Final_edgelist_cyto_low_cyto.txt",quote = F)
# write.csv(Final_edgelist_high,"Final_edgelist_cyto_high.txt",quote = F)
# 
# dim(Final_edgelist_low_cyto)
# head(Final_edgelist_low)
# 
# checkconsistency_low_raw= read.csv("Final_edgelist_cyto_low.txt",quote = T)
# checkconsistency_low = checkconsistency_low_raw[,-1]
# 
# str(checkconsistency_low$plot_order_x)
# str(Final_edgelist_low$plot_order_x)
# identical(as.numeric(checkconsistency_low$plot_order_y),Final_edgelist_low$plot_order_y)
# 
# Results_metab_cor$change$con_1

###################################################################################################
####                                    6. Differential  "expression"                       ######
#################################################################################################
library(qqplotr)
library(broom)
#
setwd("/Users/liulihe95/Desktop/metabolomics_0417")
#install.packages("gdata")
options(stringsAsFactors = FALSE)
networkData = read.xls("metabolomics_data.xlsx",sheet = 1, method=c("tab"), header= T)
rowindex = networkData[,1]
for (i in c(1:(0.5*length(networkData[,1])))){
  rowindex[i] = paste(rowindex[i],i,sep = "_")
  rowindex[i+7] = paste(rowindex[i+7],i,sep = "_")
}
rownames(networkData) = rowindex; networkData = networkData[,-1]
networkData = data.frame(t(networkData))
datExpr_metab <- as.data.frame(t(networkData));names(datExpr_metab) = rownames(networkData); rownames(datExpr_metab) = rowindex
#
raw_p = numeric(ncol(datExpr_metab))
all_qnorm = list()
all_resid = list()
TRT = substr(rownames(datExpr_metab),1,4)
for (i in c(1:ncol(datExpr_metab))){
  #test_lm = data.frame(trt = factor(index,levels = c('SF','LF')),measure=as.numeric(datExpr_metab[i,]))
  #lmod = lm(measure~.,test_lm)
  lmod = lm(datExpr_metab[,i] ~ TRT)
  raw_p[i] = summary(lmod)$coefficients[2,4]
  qqplot = ggplot(data = data.frame(residuals(lmod)), mapping = aes(sample = residuals.lmod.)) + stat_qq_band() + stat_qq_line() + stat_qq_point()+ ggtitle(allnames_label[i])+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
  all_qnorm[[i]] <- qqplot
  resid_plot = ggplot(augment(lmod), aes(x = .fitted, y = .resid)) + geom_point()+ggtitle(allnames_label[i])
  all_resid[[i]] <- resid_plot
}
### cbeck plots for signif ones
pdf("signig_QQplot.pdf")
for (i in which(raw_p< 0.05)) {
  print(all_qnorm[[i]])
}
dev.off()
pdf("signif_resid_plot.pdf")
for (i in which(raw_p< 0.05)) {
  print(all_resid[[i]])
}
dev.off()


all_resid[which(raw_p< 0.05)]

colnames
summary(raw_p)
colnames(datExpr_metab[raw_p< 0.05])
colnames(datExpr_metab[raw_p< 0.01])

# colnames(datExpr_metab[which.min(raw_p)])
# colnames(datExpr_metab[which.min(raw_p_k)])
# plot(raw_p ~ raw_p_k)

pdf("all_QQplot.pdf")
for (i in 1:ncol(datExpr_metab)) {
  print(all_qnorm[[i]])
}
dev.off()

pdf("all_resid_plot.pdf")
for (i in 1:ncol(datExpr_metab)) {
  print(all_resid[[i]])
}
dev.off()

##########   kruskal.test   ########    ####test_df = data.frame(TRT = TRT, DATA = datExpr_metab[,1])
# dim(datExpr_metab)
raw_p_k = numeric(ncol(datExpr_metab))
TRT = substr(rownames(datExpr_metab),1,4)
#kruskal.test(TRT~DATA,test_df)
for (i in c(1:ncol(datExpr_metab))){
  test_lm = data.frame(trt = substr(rownames(datExpr_metab),1,4), measure=as.numeric(datExpr_metab[,i]))
  kruskal = kruskal.test(trt ~ measure,test_lm)
  raw_p_k[i] = kruskal$p.value
  #qqplot = ggplot(data = data.frame(residuals(lmod)), mapping = aes(sample = residuals.lmod.)) + stat_qq_band() + stat_qq_line() + stat_qq_point()+ ggtitle(allnames_label[i])+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
  #all_qnorm[[i]] <- qqplot
  #resid_plot = ggplot(augment(lmod), aes(x = .fitted, y = .resid)) + geom_point()+ggtitle(allnames_label[i])
  #all_resid[[i]] <- resid_plot
}  
colnames(datExpr_metab[raw_p_k < 0.01])
colnames(datExpr_metab[raw_p_k < 0.05])
summary(raw_p_k)
# adjust p value
# adj_p = p.adjust(raw_p,method = "fdr",n = length(raw_p))
# summary(adj_p)


##########   t.test unequal var   ########   
# dim(datExpr_metab)
# ttest=t.test(datExpr_metab[(1:table(TRT)[1]),1],datExpr_metab[((table(TRT)[1]+1):length(TRT)),1],var.equal = F)
raw_p_t = numeric(ncol(datExpr_metab))
TRT = substr(rownames(datExpr_metab),1,4)
#kruskal.test(TRT~DATA,test_df)
for (i in c(1:ncol(datExpr_metab))){
  #test_lm = data.frame(trt = substr(rownames(datExpr_metab),1,4), measure=as.numeric(datExpr_metab[,i]))
  ttest = t.test(datExpr_metab[(1:table(TRT)[1]),i],datExpr_metab[((table(TRT)[1]+1):length(TRT)),i], var.equal = F)
  raw_p_t[i] = ttest$p.value
  #qqplot = ggplot(data = data.frame(residuals(lmod)), mapping = aes(sample = residuals.lmod.)) + stat_qq_band() + stat_qq_line() + stat_qq_point()+ ggtitle(allnames_label[i])+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
  #all_qnorm[[i]] <- qqplot
  #resid_plot = ggplot(augment(lmod), aes(x = .fitted, y = .resid)) + geom_point()+ggtitle(allnames_label[i])
  #all_resid[[i]] <- resid_plot
}
colnames(datExpr_metab[raw_p_t < 0.01])
colnames(datExpr_metab[raw_p_t < 0.05])
summary(raw_p_t)
#
colnames(datExpr_metab[raw_p < 0.01])
colnames(datExpr_metab[raw_p < 0.05])
which(raw_p_t < 0.05)


#############     integrate multi omics   ################
label2 = c(27:66) # Acylcarnitines
#  "C16.2.OH"  # 3-Hydroxypalmitoylcarnitine # Hydroxyhexadecadienoylcarnitine
label4 = c(81:153)  # Phosphatidylcholines #

#  "PC.aa.C36.6" "PC.ae.C36.4" 
# Phosphatidylcholine with diacyl residue sum C36:6
# Phosphatidylcholine with acyl-alkyl residue sum C36:4
mean(datExpr_low[,"C16.2.OH"])
mean(datExpr_high[,"C16.2.OH"])

mean(datExpr_low[,"PC.aa.C36.6"])
mean(datExpr_high[,"PC.aa.C36.6"])

mean(datExpr_low[,"PC.ae.C36.4"])
mean(datExpr_high[,"PC.ae.C36.4"])

mean(datExpr_low[,"C18.2"])
mean(datExpr_high[,"C18.2"])

mean(datExpr_low[,"PC.aa.C24.0"])
mean(datExpr_high[,"PC.aa.C24.0"])

mean(datExpr_low[,"PC.aa.C40.4"])
mean(datExpr_high[,"PC.aa.C40.4"])

mean(datExpr_low[,"PC.ae.C38.5"])
mean(datExpr_high[,"PC.ae.C38.5"])

mean(datExpr_low[,"SM..OH..C22.2"])
mean(datExpr_high[,"SM..OH..C22.2"])



"C18.2”

PC.aa.C24.0   
PC.aa.C40.4 
PC.ae.C38.5   

SM..OH..C22.2
