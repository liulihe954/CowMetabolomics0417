# load sda library
#install.packages("sda")
library("sda")
rew_mat1 = Results_metab_cor$cormatr_ref
rew_mat2 = Results_metab_cor$cormatr_test
adj_mat1 = Results_metab_cor$adjmatr_ref
adj_mat2 = Results_metab_cor$adjmatr_test
for(i in 1:(ncol(rew_mat)-1)){
  for (j in c((i+1):ncol(rew_mat))){
    test1 = adj_mat1[i,j]
    test2 = adj_mat2[i,j]
    if(test1 == 0){rew_mat1[i,j] = rew_mat1[j,i] = 0}
    if(test2 == 0){rew_mat2[i,j] = rew_mat2[j,i] = 0}
  }
}  
# rew_mat1 = abs(rew_mat1)
# rew_mat2 = abs(rew_mat2)
## rewiring score

all_rewire = numeric(nrow(rew_mat1))
for (i in c(1:nrow(rew_mat1))){
  vc1 = rew_mat1[i,]
  vc2 = rew_mat2[i,]
  bindmatrix = rbind(vc1,vc2)
  index = factor(rep('A',2))
  TEST = centroids(bindmatrix,index,lambda.var=0,centered.data=TRUE,verbose = 0)
  centriod = TEST$means[,2]
  rewire1 = dist(rbind(vc1, centriod),method = "euclidean")
  rewire2 = dist(rbind(vc2, centriod),method = "euclidean")
  Dn = (rewire1+rewire2)/(length(vc1)-1)
  all_rewire[i] = Dn
}
summary(all_rewire)

rownames(networkData)[rank(-abs(all_rewire),ties.method = "min") <=5]

rownames(networkData)[which.max(abs(Results_metab_cor$change$con_change))]
rownames(networkData)[which.max(abs(Results_metab_cor$change$clst_coef_change))]
rownames(networkData)[which(Results_metab_cor$change$proportion == 0)]


plot(all_rewire ~ Results_metab_cor$change$con_change)
#plot(all_rewire ~ Results_metab_cor$change$clst_coef_change)


cor(all_rewire , Results_metab_cor$change$con_change)
cor(all_rewire , Results_metab_cor$change$clst_coef_change)
summary(lm(all_rewire ~ Results_metab_cor$change$con_change))


rank(-abs(all_rewire),ties.method = "min")
rank(-abs(Results_metab_cor$change$con_change),ties.method = "min")
rank(-abs(Results_metab_cor$change$clst_coef_change),ties.method = "min")




dev.off()
## prepare data set
data(iris) # good old iris data
X = as.matrix(iris[,1:4])
Y = iris[,5]

dim(X)
length(Y)
str(iris)

## estimate centroids and empirical pooled variances
centroids(X, Y, lambda.var=0)

## also compute group-specific variances
centroids(X, Y, var.groups=TRUE, lambda.var=0)

## use shrinkage estimator for the variances
centroids(X, Y, var.groups=TRUE)

## return centered data
xc = centroids(X, Y, centered.data=TRUE)$centered.data
apply(xc, 2, mean)

## useful, e.g., to compute the inverse pooled correlation matrix
powcor.shrink(xc, alpha=-1)




