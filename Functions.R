#########################################################################################################################
## for unweighted netowrk 04062019
## get basic stats
## report cor matrix, pvalue matrix and adj matrix 
##
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
  # get proprotion of same neighboors
  proportion = numeric(ncol(datExpr1))
  for (i in c(1:ncol(datExpr1))){
    proportion[i] = (length(intersect(which(adjmatr1[i,] == 1),which(adjmatr2[i,] == 1)))-1)/(length(adjmatr1[i,] == 1)-1)
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
    proportion = proportion,
    proportion_rank = rank(-proportion,ties.method = "min"),
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


get_NC_spearman = function(datExpr1,datExpr2,r_thres,p_thres){
  options(warn = -1)
  #calcu of matrx1
  cormatr1 <- cor(datExpr1,method = c("spearman"))
  adjmatr1 = matrix(1,ncol(datExpr1),ncol(datExpr1))
  colnames(adjmatr1) = colnames(datExpr1)
  rownames(adjmatr1) = colnames(datExpr1)
  adjmatr1[abs(cormatr1) < r_thres] = 0
  #calcu of matrx2
  cormatr2 <- cor(datExpr1,method = c("spearman"))
  adjmatr2 = matrix(1,ncol(datExpr2),ncol(datExpr2))
  colnames(adjmatr2) = colnames(datExpr2)
  rownames(adjmatr2) = colnames(datExpr2)
  adjmatr2[abs(cormatr2) < r_thres] = 0
  #use threshold and get pvalue matrix
  pvalmatr1 = matrix(0,ncol(datExpr1),ncol(datExpr1))
  pvalmatr2 = matrix(0,ncol(datExpr2),ncol(datExpr2))
  for(i in 1:(ncol(datExpr1)-1)){
    for (j in c((i+1):ncol(datExpr1))){
      r1 = cor.test(datExpr1[,i],datExpr1[,j],method = c("spearman"))
      r2 = cor.test(datExpr2[,i],datExpr2[,j],method = c("spearman"))
      pvalmatr1[i,j] = pvalmatr1[j,i] = r1$p.value
      pvalmatr2[i,j] = pvalmatr2[j,i] = r2$p.value
      if(r1$p.value >= p_thres){adjmatr1[i,j] = adjmatr1[j,i] = 0}
      if(r2$p.value >= p_thres){adjmatr2[i,j] = adjmatr2[j,i] = 0}
    }
  }
  # get proprotion of same neighboors
  proportion = numeric(ncol(datExpr1))
  for (i in c(1:ncol(datExpr1))){
    proportion[i] = length(intersect(which(adjmatr1[i,] == 1),which(adjmatr2[i,] == 1)))/length(adjmatr1[1,] == 1)
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
    proportion = proportion,
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
# get_highcor_genes = function(datExpr1,r_thres,p_thres){
  #calcu of matrx1
  cormatr1 <- cor(datExpr1)
  adjmatr1 = matrix(1,ncol(datExpr1),ncol(datExpr1))
  colnames(adjmatr1) = colnames(datExpr1)
  rownames(adjmatr1) = colnames(datExpr1)
  adjmatr1[abs(cormatr1) < r_thres] = 0
  #use threshold and get pvalue matrix
  pvalmatr1 = matrix(0,ncol(datExpr1),ncol(datExpr1))
  for(i in 1:(ncol(datExpr1)-1)){
    for (j in c((i+1):ncol(datExpr1))){
      r1 = cor.test(datExpr1[,i],datExpr1[,j])
      pvalmatr1[i,j] = pvalmatr1[j,i] = r1$p.value
      if(r1$p.value >= p_thres){adjmatr1[i,j] = adjmatr1[j,i] = 0}
    }
  }
  high_index = c(rowSums(adjmatr1)-1)>0
  reduce_expression = datExpr1[,high_index]
  return(reduce_expression)
}





#########################################################################################################################
 ## rewiring score calcualtion




#########################################################################################################################
### Check zeros __ 05062019
### check a vector if it has all zeros (possibly a column of a expression dataframe)
### can set start point and stop point
### can set the threshold, how many zeros you wnat check 
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
# checkzero_result = check_zero(datExpr_low,datExpr_high,1,length(datExpr_high),7)


#########################################################################################################################
#### select gene
### select genes of con2 con2 con_change 
### every time re-rank the table/dataframe according to the corresponding column, thus we have focus
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


#########################################################################################################################
##### plot the combinded distribution plot with ggplot\
##### first put together our data (two to one) and paste the name to it as another column
##### plot parameters will be extracted from the max(data)
Plot_compare_distrib_con = function(con1,con2,title){
  meanCon_ref = mean(con1)
  meanCon_test = mean(con2)
  ref = data.frame(
    connectivity = as.numeric(con1),
    category = rep("ref",length(con1)))
  test = data.frame(
    connectivity = as.numeric(con2),
    category = rep("test",length(con2)))
  test_combine_dataset <- do.call('rbind', list(ref,test))
  ##################################        2. plotting         ##########################
  #plot3
  #par(mfrow=c(2,1))
  ggplot(test_combine_dataset, aes(x=connectivity, fill=category)) +
    geom_histogram(binwidth=1,alpha=0.6, position="identity", aes(y = ..count..), color="black") +
    # geom_density(alpha=0.6,trim = F) +
    #xlim(0,(max(c(con1,con2))+5))+
    geom_vline(aes(xintercept=meanCon_ref), color="black", linetype="dashed", size=1) +
    geom_vline(aes(xintercept=meanCon_test), color="blue", linetype="dashed", size=1) +
    theme_gray()+
    theme(legend.position="None")+
    labs(title=title, x="Connectivity", y = "Frequency")+
    theme(axis.text.x = element_text(size = 15, family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
    theme(axis.text.y = element_text(size = 15,family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
    theme(axis.title.x = element_text(size = 15,family = "Microsoft Sans Serif",color = "black",vjust = 0.5, hjust = 0.5))+
    theme(axis.title.y = element_text(size = 15, color = "black",family = "Microsoft Sans Serif", vjust = 0.5, hjust = 0.5))+
    theme(plot.title = element_text(size = 20, family = "Microsoft Sans Serif",color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
    theme(plot.title = element_text(hjust = 0.5))
}
Plot_compare_distrib_clst = function(clst1,clst2,title){
  meanClstcoef_ref = mean(clst1)
  meanClstcoef_test = mean(clst2)
    ### clst coef --- test_combine_dataset_clstcoef
  ref_clscoef = data.frame(
    clstcoef = as.numeric(clst1),
    category = rep("ref",length(clst1)))
  test_clscoef = data.frame(
    clstcoef = as.numeric(clst2),
    category = rep("test",length(clst2)))
  test_combine_dataset_clstcoef <- do.call('rbind', list(ref_clscoef,test_clscoef))
  ##################################        2. plotting         ##########################
  ggplot(test_combine_dataset_clstcoef, aes(x=clstcoef, fill=category)) +
    geom_histogram(binwidth=.01,alpha=0.6, position="identity", aes(y = ..count..), color="black") +
    #  geom_density(alpha=0.6,trim = F) +
    xlim(-0.1,1.1)+
    #ylim(0,0.1*nrow(test_combine_dataset_clstcoef))+
    geom_vline(aes(xintercept=meanClstcoef_ref), color="black", linetype="dashed", size=1) +
    geom_vline(aes(xintercept=meanClstcoef_test), color="blue", linetype="dashed", size=1) +
    theme_gray()+
    theme(legend.position="None")+
    labs(title=title, x="Cluster Coefficient", y = "Frequency")+
    theme(axis.text.x = element_text(size = 15, family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
    theme(axis.text.y = element_text(size = 15,family = "Microsoft Sans Serif",color = "black", vjust = 0.5, hjust = 0.5))+
    theme(axis.title.x = element_text(size = 15,family = "Microsoft Sans Serif",color = "black",vjust = 0.5, hjust = 0.5))+
    theme(axis.title.y = element_text(size = 15, color = "black",family = "Microsoft Sans Serif", vjust = 0.5, hjust = 0.5))+
    theme(plot.title = element_text(size = 20, family = "Microsoft Sans Serif",color = "black", face = "bold", vjust = 0.5, hjust = 0.5))+
    theme(plot.title = element_text(hjust = 0.5))
}


#########################################################################################################################
### plot multiple figures in one
###
###cols for columns
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#########################################################################################################################
### sometime GetNetworkConcept does not work (dont konw why, so weird) /// in WGCNA 
### look into the source code and find the functionit self
### tear them apart for further use
get_clst_coef = function (adjmat1) {
  no.nodes = dim(adjmat1)[[1]]
  computeLinksInNeighbors <- function(x, imatrix) {
    x %*% imatrix %*% x
  }
  computeSqDiagSum = function(x, vec) {
    sum(x^2 * vec)
  }
  nolinksNeighbors <- c(rep(-666, no.nodes))
  total.edge <- c(rep(-666, no.nodes))
  maxh1 = max(as.dist(adjmat1))
  minh1 = min(as.dist(adjmat1))
  if (maxh1 > 1 | minh1 < 0) {
    stop(paste("ERROR: the adjacency matrix contains entries that are larger", 
               "than 1 or smaller than 0: max=", maxh1, ", min=", 
               minh1))
  }
  else {
    nolinksNeighbors <- apply(adjmat1, 1, computeLinksInNeighbors, 
                              imatrix = adjmat1)
    subTerm = apply(adjmat1, 1, computeSqDiagSum, vec = diag(adjmat1))
    plainsum <- apply(adjmat1, 1, sum)
    squaresum <- apply(adjmat1^2, 1, sum)
    total.edge = plainsum^2 - squaresum
    CChelp = rep(-666, no.nodes)
    CChelp = ifelse(total.edge == 0, 0, (nolinksNeighbors - subTerm)/total.edge)
    CChelp
  }
}

get_DCH = function (adj) { Size = dim(adj)[1]
  Connectivity = apply(adj, 2, sum)
  Density = sum(Connectivity)/(Size * (Size - 1))
  Clst_coef = get_clst_coef(adj)
  Centralization = Size * (max(Connectivity) - mean(Connectivity))/((Size - 1) * (Size - 2))
  Heterogeneity = sqrt(Size * sum(Connectivity^2)/sum(Connectivity)^2 -  1)
  Results = list(Connectivity = Connectivity, Density = Density, ClusterCoef = Clst_coef , Centralization = Centralization, Heterogeneity = Heterogeneity)
  return(Results)
}

# test_DCH
get_NC_cor_modf = function(datExpr1,datExpr2,r_thres,p_thres){
  #calcu of matrx1
  cormatr1 <- cor(datExpr1)
  adjmatr1 = matrix(1,ncol(datExpr1),ncol(datExpr1))
  colnames(adjmatr1) = colnames(datExpr1);rownames(adjmatr1) = colnames(datExpr1)
  adjmatr1[abs(cormatr1) < r_thres] = 0
  #calcu of matrx2
  cormatr2 <- cor(datExpr2)
  adjmatr2 = matrix(1,ncol(datExpr2),ncol(datExpr2))
  colnames(adjmatr2) = colnames(datExpr2);rownames(adjmatr2) = colnames(datExpr2)
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
  NC1 = get_DCH(adjmatr1)
  NC2 = get_DCH(adjmatr2)
  #combine, rank and show
  basic_results = data.frame(density1 = NC1$Density,
                             density2 = NC2$Density,
                             centralization1 = NC1$Centralization,
                             centralization2 = NC2$Centralization,
                             heterogeneity1 = NC1$Heterogeneity,
                             heterogeneity2 = NC2$Heterogeneity)
  change_results = data.frame(#gene = colnames(datExpr1),
    con_1 = NC1$Connectivity,
    #scl_con_1 = (NC1$Connectivity)/max(NC1$Connectivity),
    con_2 = NC2$Connectivity,
    #scl_con_2 = (NC2$Connectivity)/max(NC2$Connectivity),
    con_change = -(NC1$Connectivity - NC2$Connectivity),
    #scl_con_change = -(scl_con_1  - scl_con_2),
    rank_scl_con = rank(-abs(NC1$Connectivity-NC2$Connectivity),ties.method = "min"),
    cls_coef_1 = NC1$ClusterCoef,
    cls_coef_2 = NC2$ClusterCoef,
    clst_coef_change = c(NC1$ClusterCoef - NC2$ClusterCoef),
    rank_clstcoef = rank(-abs(NC1$ClusterCoef-NC2$ClusterCoef),ties.method = "min"))
  Results = list(NC1 = NC1, NC2 = NC2, cormatr_ref = cormatr1, cormatr_test = cormatr2, pvalue1 = pvalmatr1,pvalue2 = pvalmatr2 ,adjmatr_ref = adjmatr1, adjmatr_test = adjmatr2, basic = basic_results,change = change_results)
  return(Results)
}

get_NC_spearman_modf = function(datExpr1,datExpr2,r_thres,p_thres){
  options(warn = -1)
  #calcu of matrx1
  cormatr1 <- cor(datExpr1,method = c("spearman"))
  adjmatr1 = matrix(1,ncol(datExpr1),ncol(datExpr1))
  colnames(adjmatr1) = colnames(datExpr1)
  rownames(adjmatr1) = colnames(datExpr1)
  adjmatr1[abs(cormatr1) < r_thres] = 0
  #calcu of matrx2
  cormatr2 <- cor(datExpr1,method = c("spearman"))
  adjmatr2 = matrix(1,ncol(datExpr2),ncol(datExpr2))
  colnames(adjmatr2) = colnames(datExpr2)
  rownames(adjmatr2) = colnames(datExpr2)
  adjmatr2[abs(cormatr2) < r_thres] = 0
  #use threshold and get pvalue matrix
  pvalmatr1 = matrix(0,ncol(datExpr1),ncol(datExpr1))
  pvalmatr2 = matrix(0,ncol(datExpr2),ncol(datExpr2))
  for(i in 1:(ncol(datExpr1)-1)){
    for (j in c((i+1):ncol(datExpr1))){
      r1 = cor.test(datExpr1[,i],datExpr1[,j],method = c("spearman"))
      r2 = cor.test(datExpr2[,i],datExpr2[,j],method = c("spearman"))
      pvalmatr1[i,j] = pvalmatr1[j,i] = r1$p.value
      pvalmatr2[i,j] = pvalmatr2[j,i] = r2$p.value
      if(r1$p.value >= p_thres){adjmatr1[i,j] = adjmatr1[j,i] = 0}
      if(r2$p.value >= p_thres){adjmatr2[i,j] = adjmatr2[j,i] = 0}
    }
  }
  #get all the basic NC
  NC1 = get_DCH(adjmatr1)
  NC2 = get_DCH(adjmatr2)
  #combine, rank and show
  basic_results = data.frame(density1 = NC1$Density,
                             density2 = NC2$Density,
                             centralization1 = NC1$Centralization,
                             centralization2 = NC2$Centralization,
                             heterogeneity1 = NC1$Heterogeneity,
                             heterogeneity2 = NC2$Heterogeneity)
  change_results = data.frame(#gene = colnames(datExpr1),
    con_1 = NC1$Connectivity,
    #scl_con_1 = (NC1$Connectivity)/max(NC1$Connectivity),
    con_2 = NC2$Connectivity,
    #scl_con_2 = (NC2$Connectivity)/max(NC2$Connectivity),
    con_change = -(NC1$Connectivity - NC2$Connectivity),
    #scl_con_change = -(scl_con_1  - scl_con_2),
    rank_scl_con = rank(-abs(NC1$Connectivity-NC2$Connectivity),ties.method = "min"),
    cls_coef_1 = NC1$ClusterCoef,
    cls_coef_2 = NC2$ClusterCoef,
    clst_coef_change = c(NC1$ClusterCoef - NC2$ClusterCoef),
    rank_clstcoef = rank(-abs(NC1$ClusterCoef-NC2$ClusterCoef),ties.method = "min"))
  Results = list(NC1 = NC1, NC2 = NC2, cormatr_ref = cormatr1, cormatr_test = cormatr2, pvalue1 = pvalmatr1,pvalue2 = pvalmatr2 ,adjmatr_ref = adjmatr1, adjmatr_test = adjmatr2, basic = basic_results,change = change_results)
  return(Results)
}

#########################################################################################################################
### ggplot
### to plot qqplot
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}


###qq plot in ggplot: input must be data.frame/// e.g. # smp <- data.frame(norm = rnorm(100))
### package: install.packages("qqplotr") ; library(qqplotr)
# qqplot = ggplot(data = data.frame(residuals(lmod)), mapping = aes(sample = residuals.lmod.)) + stat_qq_band() + stat_qq_line() + stat_qq_point()+ ggtitle(allnames_label[i])+labs(x = "Theoretical Quantiles", y = "Sample Quantiles")


### install.packages("car");library(car)
### qqplot 
#e.g.
#
#  qqPlot(residuals(lmod))
#