##### library #####
# library(WGCNA);library(ppcor);library(igraph);library(gdata);library(ggplot2)
# library(ggpubr);require(cowplot);library(extrafont);library(dplyr);library(plotly) ;library(geomnet);library(readxl)
##################################################################################################
###                                       1. dataprep                                      ######
################################################################################################
#install.packages("gdata")
options(stringsAsFactors = FALSE)
networkData = read.xls("metabolomics_data.xlsx",sheet = 1, method=c("tab"), header= TRUE)

rowindex = networkData[,1]
for (i in c(1:(0.5*length(networkData[,1])))){
  rowindex[i] = paste(rowindex[i],i,sep = "_")
  rowindex[i+7] = paste(rowindex[i+7],i,sep = "_")
}
rownames(networkData) = rowindex; networkData = networkData[,-1]

networkData = data.frame(t(networkData))
datExpr_metab <- as.data.frame(t(networkData));names(datExpr_metab) = rownames(networkData); rownames(datExpr_metab) = rowindex
#dim(datExpr_metab)
#table(colnames(networkData) == rownames(datExpr_metab))
datExpr_low <- datExpr_metab[c(1:7),];datExpr_high <- datExpr_metab[c(8:14),]
# dim(datExpr_low);dim(datExpr_high)
###########
### making labels for sub-net --- different chemical properties
allnames_label = colnames(datExpr_low)
label1 = c(1:26) # Amino acids and Biogenic Amines 
label2 = c(27:66) # Acylcarnitines 
label3 = c(1:80)  # Lysophosphatidylcholines
label4 = c(1:153)  # Phosphatidylcholines
label5 = c(1:167)  # Sphingomyelins 
label6 = c(168)    # Hexoses 
label7 = c(167:179) # Eicosanoids & oxidation products of polyunsaturated fatty acids
#########
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
# Results_metab_cor 

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
topgene_metab_con = data.frame(
  Con_ref = Results_metab_cor$change$con_1,
  Con_rank_ref = rank(-Results_metab_cor$change$con_1,ties.method = "min"),
  #Con_ref_scl = Results_sperm_cor$change$scl_con_1,
  Con_test = Results_metab_cor$change$con_2,
  #Con_ref_scl = Results_sperm_cor$change$scl_con_1,
  Con_rank_test = rank(-Results_metab_cor$change$con_2,ties.method = "min"),
  Con_Change = Results_metab_cor$change$con_change,
  ConChange_rank =rank(-(abs(Results_metab_cor$change$con_1 - Results_metab_cor$change$con_2)),ties.method = "min")
)
rownames(topgene_metab_con) = colnames(datExpr_high_new)
######### assemble dataset 2 _ clscoef
topgene_metab_clscoef = data.frame(
  Clscoef_ref = Results_metab_cor$change$cls_coef_1 ,
  Clscoef_rank_ref = rank(-Results_metab_cor$change$cls_coef_1,ties.method = "min"),
  Clscoef_test = Results_metab_cor$change$cls_coef_2,
  Clscoef_rank_test = rank(-Results_metab_cor$change$cls_coef_2,ties.method = "min"),
  Clscoef_Change = Results_metab_cor$change$clst_coef_change,
  ConChange_rank =  rank(-(abs(Results_metab_cor$change$cls_coef_1 - Results_metab_cor$change$cls_coef_2)),ties.method = "min")
)
rownames(topgene_metab_clscoef) = colnames(datExpr_high_new)

#### define function 
SelectGene_un_cor = function(dataset,topnumber){
  index1 = dataset[,2] %in% c(1:topnumber)
  index2 = dataset[,4] %in% c(1:topnumber)
  index3 = dataset[,6] %in% c(1:topnumber)
  summary1 = dataset[index1,];summary1 = summary1[order(summary1[,2]),]
  summary2 = dataset[index2,];summary2 = summary2[order(summary2[,4]),]
  summary3 = dataset[index3,];summary3 = summary3[order(summary3[,6]),]
  summary = list(
    ref = summary1,
    test = summary2,
    change =summary3
  )
  return(summary)
}
##### example result, both in top 22 - con and clustercoef
top10_con_change = rownames(SelectGene_un_cor(topgene_metab_con,22)[[3]])
top10_ClusterCoef_change = rownames(SelectGene_un_cor(topgene_metab_clscoef,22)[[3]])
intersect(top10_con_change,top10_ClusterCoef_change)
Results_metab_cor$basic

###################################################################################################
####                                       3. Plotting                                      ######
#################################################################################################
#################################        1. generating dataset ####################
#####  con --- test_combine_dataset
ref = data.frame(
  connectivity = as.numeric(Results_metab_cor$NC1$fundamentalNCs$Connectivity),
  category = rep("ref",ncol(datExpr_low_new)))
test = data.frame(
  connectivity = as.numeric(Results_metab_cor$NC2$fundamentalNCs$Connectivity),
  category = rep("test",ncol(datExpr_high_new)))
test_combine_dataset <- do.call('rbind', list(ref,test))
str(test_combine_dataset)
# table(test_combine_dataset$category)

### clst coef --- test_combine_dataset_clstcoef
ref_clscoef = data.frame(
  clstcoef = as.numeric(Results_metab_cor$NC1$fundamentalNCs$ClusterCoef),
  category = rep("ref",ncol(datExpr_low_new)))
test_clscoef = data.frame(
  clstcoef = as.numeric(Results_metab_cor$NC2$fundamentalNCs$ClusterCoef),
  category = rep("test",ncol(datExpr_high_new)))
test_combine_dataset_clstcoef <- do.call('rbind', list(ref_clscoef,test_clscoef))
str(test_combine_dataset_clstcoef)

##################################        2. plotting         ##########################

#plot3
ggplot(test_combine_dataset, aes(x=connectivity, fill=category)) +
  geom_histogram(binwidth=1,alpha=0.6, position="identity", aes(y = ..count..), color="black") +
  # geom_density(alpha=0.6,trim = F) +
  xlim(0,80)+
  geom_vline(aes(xintercept=meanCon_ref), color="black", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=meanCon_test), color="blue", linetype="dashed", size=1) +
  theme_gray()+
  theme(legend.position="None")+
  labs(title="Distribution of Connectivity", x="Connectivity", y = "Frequency")+
  theme(axis.text.x = element_text(size = 15, family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 15,family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x = element_text(size = 15,family = "Microsoft Sans Serif",color = "black",vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 15, color = "black",family = "Microsoft Sans Serif", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(size = 20, family = "Microsoft Sans Serif",color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))

#  plot4 
ggplot(test_combine_dataset_clstcoef, aes(x=clstcoef, fill=category)) +
  geom_histogram(binwidth=.01,alpha=0.6, position="identity", aes(y = ..count..), color="black") +
  #  geom_density(alpha=0.6,trim = F) +
 # xlim(-0.01,1.2)+
  geom_vline(aes(xintercept=meanClstcoef_ref), color="black", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=meanClstcoef_test), color="blue", linetype="dashed", size=1) +
  theme_gray()+
  theme(legend.position="None")+
  labs(title="Distribution of Cluster Coefficient", x="Cluster Coefficient", y = "Frequency")+
  theme(axis.text.x = element_text(size = 15, family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 15,family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.x = element_text(size = 15,family = "Microsoft Sans Serif",color = "black",vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 15, color = "black",family = "Microsoft Sans Serif", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(size = 20, family = "Microsoft Sans Serif",color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5))

#tiff("Figure1_DistributionChange3.tiff", width = 14, height = 12, units = 'in', res = 300)
#plot_grid(plot3, plot4, align = c("v"),labels = c("A","B"), nrow = 2,label_size= 20, label_colour = "darkgreen")
dev.off()



#### testing 
dice_results <- c(1,3,2,4,5,6,5,3,2,1,6,2,6,5,6,4,1,3,2,4,6,4,1,6,3,2,4,3,4,5,6,7,1)
ggplot(dice_results,aes(x =(dice_results))) + geom_histogram(binwidth=1, colour="black", fill="white")
