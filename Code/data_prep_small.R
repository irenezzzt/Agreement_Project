
###################################
####    Agreement Project    ######
#### Rcode 1 Data Preparation #####
###################################


library(FactoMineR)
library(factoextra)
library(xtable)
library(corrplot)

# read in data and pre-process ---------------
#PATH <- "~/SCRATCH/khaoula.belahsen"

#load(file.path(PATH, "AGREEMENT/data/MH.SARDs.Novoalign.Modeling.All_Probes_Data.V2.RData"))


setwd("~/Desktop/Agreement Project/")
load(file.path("data_prep/MH.SARDs.Novoalign.Modeling.All_Probes_Data.V2.RData"))

# write.table(colnames(training),"colname.txt")
# colnames(training)

colnames(training)[18]=colnames(testing)[18]="colour_channel"
colnames(training)[21]=colnames(testing)[21]="island_relation"
rownames(training)=seq(1,nrow(training))
rownames(testing)=seq(1,nrow(testing))

set.seed(123)
# remove rows with NA values.
keep=apply(training,1,function(x) all(!is.na(x))) 
training=training[keep,]
keep=apply(testing,1,function(x) all(!is.na(x)))
testing=testing[keep,]
rm(keep)

# testing=as.data.frame(testing)
# load(file.path(PATH,"AGREEMENT/models/RF_100.RData"))
# load(file.path(PATH,"AGREEMENT/models/xgb100.RData"))
# load(file.path(PATH,"AGREEMENT/data/testagg.RData"))
# load(file.path(PATH,"AGREEMENT/data/predictions.RData"))


# remove indiv with u+m > 100
sub.training<- training[ (training[,"meth_count"]+training[,"unmeth_count"])<=100, ]
# dim(sub.training)[1]-dim(training)[1] 19949 rows removed

#training[,11:16] : CG.F60    CG.R60  CG.F100 CG.R100   CG.F120   CG.R120, CG content
count_m <- cbind(sub.training[,2],sub.training[,"meth_count"],1,sub.training[,11:16])
count_u <- cbind(sub.training[,2],sub.training[,"unmeth_count"],0,sub.training[,11:16])

colnames(count_m)<-colnames(count_u)<-c("id","count","meth","CG.F60","CG.R60","CG.F100","CG.R100","CG.F120","CG.R120")


# take the first 5000 rows (1 individuals)
set.seed(1227)
rd=sample(dim(sub.training)[1], 5000, replace = F)


train <- rbind(count_m[rd,],count_u[rd,])
train = as.data.frame(train)




# this part has problem running 

#testing=as.data.frame(testing)
#testing$num<- 1:nrow(testing) 
#testing=testing[testing$num %in% test_m$num,]


# removing rows with u+m >100
sub.testing<-testing[(testing[,"meth_count"]+testing[,"unmeth_count"])<=100,]
# dim(sub.testing)-dim(testing) 1430 rows removed

count_m <- cbind(sub.testing[,2],sub.testing[,"meth_count"],1,sub.testing[,11:16])
count_u <- cbind(sub.testing[,2],sub.testing[,"unmeth_count"],0,sub.testing[,11:16])
colnames(count_m) <-colnames(count_u) <- c("id","count","meth","CG.F60","CG.R60","CG.F100","CG.R100","CG.F120","CG.R120")
set.seed(1005)
rd<-sample( dim(count_m)[1], 5000, replace = F)
test_m=as.data.frame(count_m[rd,])
test_u=as.data.frame(count_u[rd,])

rm(count_m); rm(count_u);rm(rd)

# check range of counts in training and test set
summary(train$count) 
# heavily skewed, 3rd quater 28, max 322
summary(test_u$count); summary(test_m$count) 


train$num=1:nrow(train)
test_m$num=1:nrow(test_m)
test_u$num=1:nrow(test_u)

#test_m=test_m[test_m$count < 100,]
#test_u=test_u[test_u$num %in% test_m$num,]
#test_u=test_u[test_u$count < 100,]
#test_m=test_m[test_m$num %in% test_u$num,]


#save(train,test_m,test_u,file=file.path(PATH,"AGREEMENT/data/testagg.RData"))
save(train,test_m,test_u,file=file.path("data_prep/testagg.RData"))

# training PCA model on the GC content -------------

#correlation?
corrplot(cor(train[,c(2,4:9)]), method="circle")
# colnames(train[,3:8]) #F60-R120

pca.tr <-PCA(train[,4:9])
eig.val <- pca.tr$eig # first 2 gives 96% of the variance explained

xtable(eig.val, digits=4)
fviz_eig(pca.tr, addlabels = TRUE, ylim = c(0, 60))
fviz_pca_var(pca.tr, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Évite le chevauchement de texte
)



# I don't think we should do PCA on test set. ------
# We should use the projection matrix from training set PCA and extract the same PCs from test set ---------
# Hence No Need to Run the Following: ---------
# PCA test_m 

res.pca.tm <-PCA(test_m[,3:8])
eig.val <- res.pca.tm$eig

fviz_eig(res.pca.tm, addlabels = TRUE, ylim = c(0, 60))
fviz_pca_var(res.pca.tm, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Évite le chevauchement de texte
)

# Contributions des variables à PC1
fviz_contrib(res.pca.tm, choice = "var", axes = 1, top = 10)

# Contributions des variables à PC2
fviz_contrib(res.pca.tm, choice = "var", axes = 2, top = 10)

# Contributions à PC1 et PC2 
fviz_contrib(res.pca.tm, choice = "var", axes = 1:2, top = 10)

#PC1 and PC2 

test_m$PC1=res.pca.tm$ind$coord[,1]
test_m$PC2=res.pca.tm$ind$coord[,2]

# PCA test_u

res.pca.tu <-PCA(test_u[,3:8])
eig.val <- res.pca.tu$eig
fviz_eig(res.pca.tu, addlabels = TRUE, ylim = c(0, 60))
fviz_pca_var(res.pca.tu, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Évite le chevauchement de texte
)

# Contributions des variables à PC1
fviz_contrib(res.pca.tu, choice = "var", axes = 1, top = 10)

# Contributions des variables à PC2
fviz_contrib(res.pca.tu, choice = "var", axes = 2, top = 10)

# Contributions à PC1 et PC2 
fviz_contrib(res.pca.tu, choice = "var", axes = 1:2, top = 10)

#PC1 and PC2 

test_u$PC1=res.pca.tu$ind$coord[,1]
test_u$PC2=res.pca.tu$ind$coord[,2]

save(pca.tr,res.pca.tm,res.pca.tu,file=file.path("data_prep/modpca.RData"))
save(train,test_m,test_u,file=file.path("data_prep/datapca.RData"))
save(training,testing,file=file.path("data_prep/dataing.RData"))    









# probaly check if PCs are similar? -------
res.pca.tm <-prcomp(test_m[,4:9])
self_pcs<-predict(res.pca.tm, test_m[,4:9])
pred_pcs<-predict(pca.tr, test_m[,4:9])

head(pred_pcs[[1]])
diag(cor(self_pcs, pred_pcs))
corrplot(cor(self_pcs, pred_pcs))
test_m$PC1=res.pca.tm$ind$coord[,1]
test_m$PC2=res.pca.tm$ind$coord[,2]

# save some new forms of data set --------
# save PC1, PC2, and their polynomials (why?)

pca.tr<-prcomp(train[,4:9])

train$PC1=predict(pca.tr, test_m[,4:9])[,1]
train$PC2=predict(pca.tr, test_m[,4:9])[,2]
train$PC12=(train$PC1)^2
train$PC22=(train$PC2)^2
train$PC13=(train$PC1)^3
train$PC23=(train$PC2)^3

test.pcs<-predict(pca.tr, test_m[,4:9])

test_m$PC1 <- test.pcs[,1]
test_m$PC2 <- test.pcs[,2]
test_m$PC12=(test_m$PC1)^2
test_m$PC22=(test_m$PC2)^2
test_m$PC13=(test_m$PC1)^3
test_m$PC23=(test_m$PC2)^3

test.pcs<-predict(pca.tr, test_u[,4:9])

test_u$PC1 <- test.pcs[,1]
test_u$PC2 <- test.pcs[,2]
test_u$PC12=(test_u$PC1)^2
test_u$PC22=(test_u$PC2)^2
test_u$PC13=(test_u$PC1)^3
test_u$PC23=(test_u$PC2)^3

# "s" set , what does s mean?
train_s=train
testm_s=test_m
testu_s=test_u

train_s$CG.F602=(train_s$CG.F60)^2
train_s$CG.F603=(train_s$CG.F60)^3
train_s$CG.R602=(train_s$CG.R60)^2
train_s$CG.R603=(train_s$CG.R60)^3
train_s$CG.F1002=(train_s$CG.F100)^2
train_s$CG.F1003=(train_s$CG.F100)^3
train_s$CG.R1002=(train_s$CG.R100)^2
train_s$CG.R1003=(train_s$CG.R100)^3
train_s$CG.F1202=(train_s$CG.F120)^2
train_s$CG.F1203=(train_s$CG.F120)^3
train_s$CG.R1202=(train_s$CG.R120)^2
train_s$CG.R1203=(train_s$CG.R120)^3
testm_s$CG.F602=(testm_s$CG.F60)^2
testm_s$CG.F603=(testm_s$CG.F60)^3
testm_s$CG.R602=(testm_s$CG.R60)^2
testm_s$CG.R603=(testm_s$CG.R60)^3
testm_s$CG.F1002=(testm_s$CG.F100)^2
testm_s$CG.F1003=(testm_s$CG.F100)^3
testm_s$CG.R1002=(testm_s$CG.R100)^2
testm_s$CG.R1003=(testm_s$CG.R100)^3
testm_s$CG.F1202=(testm_s$CG.F120)^2
testm_s$CG.F1203=(testm_s$CG.F120)^3
testm_s$CG.R1202=(testm_s$CG.R120)^2
testm_s$CG.R1203=(testm_s$CG.R120)^3


testu_s$CG.F602=(testu_s$CG.F60)^2
testu_s$CG.F603=(testu_s$CG.F60)^3
testu_s$CG.R602=(testu_s$CG.R60)^2
testu_s$CG.R603=(testu_s$CG.R60)^3
testu_s$CG.F1002=(testu_s$CG.F100)^2
testu_s$CG.F1003=(testu_s$CG.F100)^3
testu_s$CG.R1002=(testu_s$CG.R100)^2
testu_s$CG.R1003=(testu_s$CG.R100)^3
testu_s$CG.F1202=(testu_s$CG.F120)^2
testu_s$CG.F1203=(testu_s$CG.F120)^3
testu_s$CG.R1202=(testu_s$CG.R120)^2
testu_s$CG.R1203=(testu_s$CG.R120)^3

# factorize the methylated counts is not realistic thus is not considered

save(train_s,testm_s,testu_s,file=file.path("data_prep/data_s.RData"))





