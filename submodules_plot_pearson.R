label3 = c(67:80)
datExpr_low_sub3 = datExpr_low[,label3]; datExpr_high_sub3 = datExpr_high[,label3]

### only for subset 3
  #calcu of matrx1
  cormatr1 <- cor(datExpr_low_sub3)
  adjmatr1 = matrix(1,ncol(datExpr_low_sub3),ncol(datExpr_low_sub3))
  colnames(adjmatr1) = colnames(datExpr_low_sub3)
  rownames(adjmatr1) = colnames(datExpr_low_sub3)
  adjmatr1[abs(cormatr1) < 0.5] = 0
  #calcu of matrx2
  cormatr2 <- cor(datExpr_high_sub3)
  adjmatr2 = matrix(1,ncol(datExpr_high_sub3),ncol(datExpr_high_sub3))
  colnames(adjmatr2) = colnames(datExpr_high_sub3)
  rownames(adjmatr2) = colnames(datExpr_high_sub3)
  adjmatr2[abs(cormatr2) < 0.5] = 0
  #use threshold and get pvalue matrix
  pvalmatr1 = matrix(0,ncol(datExpr_low_sub3),ncol(datExpr_low_sub3))
  pvalmatr2 = matrix(0,ncol(datExpr_high_sub3),ncol(datExpr_high_sub3))
  for(i in 1:(ncol(datExpr_low_sub3)-1)){
    for (j in c((i+1):ncol(datExpr_low_sub3))){
      r1 = cor.test(datExpr_low_sub3[,i],datExpr_low_sub3[,j])
      r2 = cor.test(datExpr_high_sub3[,i],datExpr_high_sub3[,j])
      pvalmatr1[i,j] = pvalmatr1[j,i] = r1$p.value
      pvalmatr2[i,j] = pvalmatr2[j,i] = r2$p.value
      if(r1$p.value >= 0.05){adjmatr1[i,j] = adjmatr1[j,i] = 0}
      if(r2$p.value >= 0.05){adjmatr2[i,j] = adjmatr2[j,i] = 0}
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

  
getAnywhere(.ClusterCoef.fun)
  
test_cls = get_clst_coef(adjmatr2)
??.ClusterCoef.fun()
 

