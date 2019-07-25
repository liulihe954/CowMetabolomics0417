##### library #####
library(WGCNA);library(ppcor);library(gdata)
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
#dim(datExpr_metab)
#table(colnames(networkData) == rownames(datExpr_metab))
datExpr_low <- datExpr_metab[c(1:7),];datExpr_high <- datExpr_metab[c(8:14),]
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
label3 = c(67:80)
test_data_sub3_low = datExpr_low[,label3]; test_data_sub3_high = datExpr_high[,label3]
testResults_cor_sub3 = get_NC_cor(test_data_sub3_low,test_data_sub3_high,0.5,0.05)
save(testResults_cor_sub3, file = "testResults_cor_sub3.RData")
