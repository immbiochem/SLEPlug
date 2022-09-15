
library(ggplot2)
library(dplyr)
library(xlsx)
library(lattice)
library(vegan)
library(mixOmics)
set.seed(30) 
library(MASS)
library(MLmetrics)
library(tree)
library(randomForest)
library(pROC)
library(car)


# Download the data
SLE = read.xlsx("SLE_metabolites.xlsx", sheetIndex = 1)
str(SLE)
#


# NDEF ~ 0
for (i in c(1:ncol(SLE))){
  SLE[,i] = sapply(SLE[,i], function(x) ifelse(x=='NDEF', 0, x))
}

# TAG~NA
for (i in c(1:ncol(SLE))){
  SLE[,i] = sapply(SLE[,i], function(x) ifelse(x=='TAG', NA, x))
}

# make it numeric
make_it_num = function(df, q){
  for (i in q) {
    df[,i] = as.numeric(df[,i])
  }
  return (df)
}
SLE = make_it_num(SLE, c(3:230)) 
#
SLE$Phenotype = factor(SLE$Phenotype)
histogram(SLE$Phenotype)

# NA imput
#
na_counts = sapply(SLE[1:39,], 
                   function(x) sum(sapply(x, function(z) is.na(z)))/length(1:39))
na_counts[which(na_counts>0.1)]
# none
#
na_counts_sle_np = sapply(SLE[40:89,], 
                   function(x) sum(sapply(x, function(z) is.na(z)))/length(40:89))
na_counts_sle_np[which(na_counts_sle_np>0.1)]
#
# Pyr 0.26
na_counts_sle_p = sapply(SLE[90:119,], 
                          function(x) sum(sapply(x, function(z) is.na(z)))/length(90:119))
na_counts_sle_p[which(na_counts_sle_p>0.1)]
#
# Pyr 0.23
#
# Pyr has a lot of NA's in SLE patients, so we remove it
SLE = SLE[,-which(colnames(SLE) == 'Pyr')]
#

# Use median-imput
minor_na = which(sapply(SLE, function(x) sum(is.na(x)))>0)
#
for (i in minor_na){
  hc_median = median(SLE[,i][which(SLE$Phenotype == 'HC')], na.rm=T)
  sle_np_median = median(SLE[,i][which(SLE$Phenotype == 'SLE-NP')], na.rm=T)
  sle_p_median = median(SLE[,i][which(SLE$Phenotype == 'SLE-P')], na.rm=T)
  SLE[,i] = c(sapply(SLE[,i][1:39], function(x) ifelse(is.na(x), hc_median, x)),
                 sapply(SLE[,i][40:89], function(x) ifelse(is.na(x), sle_np_median, x)),
              sapply(SLE[,i][90:119], function(x) ifelse(is.na(x), sle_p_median, x)))
}
# Chech inputation
sum(is.na(SLE)) # 0

#
#write.xlsx(SLE, "SLE_transfor.xlsx")

# _____________________________________________________________________________
#
# ALL groups
ggplot(data.frame(SLE[,c(2:229)], scores(rda(SLE[,c(3:229)], scale = TRUE),
                                                  display = "sites", choices = c(1, 2, 3), 
                                                  scaling = "sites")), 
       aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Phenotype), alpha = 0.8) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle(label = "Ординация в осях главных компонент") + 
  theme_bw()

# HC vs SLE-NP
ggplot(data.frame(SLE[1:89,c(2:229)], scores(rda(SLE[1:89,c(3:229)], scale = TRUE),
                                         display = "sites", choices = c(1, 2, 3), 
                                         scaling = "sites")), 
       aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Phenotype), alpha = 0.8) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle(label = "Ординация в осях главных компонент") + 
  theme_bw()


# HC vs SLE-P
ggplot(data.frame(SLE[c(1:39, 90:119),c(2:229)], scores(rda(SLE[c(1:39, 90:119),c(3:229)], scale = TRUE),
                                             display = "sites", choices = c(1, 2, 3), 
                                             scaling = "sites")), 
       aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Phenotype), alpha = 0.8) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle(label = "Ординация в осях главных компонент") + 
  theme_bw()

# SLE-P vs SLE-NP
ggplot(data.frame(SLE[40:119,c(2:229)], scores(rda(SLE[40:119,c(3:229)], scale = TRUE),
                                             display = "sites", choices = c(1, 2, 3), 
                                             scaling = "sites")), 
       aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Phenotype), alpha = 0.8) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle(label = "Ординация в осях главных компонент") + 
  theme_bw()


### Take into account amino acids only

#  _________________________________________________________________________
# REDUCTION
cap_met = c('AcAce','Ace','Ala','Alb',
  'bOHBut','Cit','Crea',
  'Glc','Gln','Glol','Gly',
  'His','Ile','Lac',
  'Leu','Phe','Tyr','Val')

SLE_red = SLE[c('Phenotype',cap_met)]

#
# ALL groups
ggplot(data.frame(SLE_red, scores(rda(SLE_red[,c(2:19)], scale = TRUE),
                                         display = "sites", choices = c(1, 2, 3), 
                                         scaling = "sites")), 
       aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Phenotype), alpha = 0.8) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle(label = "Ординация в осях главных компонент") + 
  theme_bw()

# HC vs SLE-NP
ggplot(data.frame(SLE_red[1:89,], scores(rda(SLE_red[1:89,c(3:19)], scale = TRUE),
                                             display = "sites", choices = c(1, 2, 3), 
                                             scaling = "sites")), 
       aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Phenotype), alpha = 0.8) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle(label = "Ординация в осях главных компонент") + 
  theme_bw()


# HC vs SLE-P
ggplot(data.frame(SLE_red[c(1:39, 90:119),], scores(rda(SLE_red[c(1:39, 90:119),c(3:19)], scale = TRUE),
                                                        display = "sites", choices = c(1, 2, 3), 
                                                        scaling = "sites")), 
       aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Phenotype), alpha = 0.8) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle(label = "Ординация в осях главных компонент") + 
  theme_bw()

# SLE-P vs SLE-NP
ggplot(data.frame(SLE_red[40:119,], scores(rda(SLE_red[40:119,c(3:19)], scale = TRUE),
                                               display = "sites", choices = c(1, 2, 3), 
                                               scaling = "sites")), 
       aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Phenotype), alpha = 0.8) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) + 
  ggtitle(label = "Ординация в осях главных компонент") + 
  theme_bw()


#____________________________________________________________________________________

# PLS-DA

# ALL DATA ALL GROUPS
list.keepX <- c(5:10,  seq(20, 100, 10))

tune.splsda.srbct <- tune.splsda(SLE[,c(3:229)], SLE$Phenotype, ncomp = 2, # число компонент берется как число групп-1
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 10)   
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp 
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp] 
#Оптимизированный sPLS-DA
MyResult.splsda <- splsda(SLE[,c(3:229)], SLE$Phenotype, 
                          ncomp = 2, keepX = select.keepX)
#
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(MyResult.splsda, comp = 1, 
             size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", 
             contrib = "max", legend = FALSE, col.ties="black", ndisplay = 5)
plotLoadings(MyResult.splsda, comp = 2, 
             size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", 
             contrib = "max",ndisplay = 5,  legend = FALSE, col.ties="black")
plotIndiv(MyResult.splsda, ind.names = F, 
          ellipse = T, style = "graphics", abline = TRUE, cex = 2, 
          pch = 19, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, 
          title = "sPLS-DA ordination of samples", size.title = 1.5)
legend("bottomright", legend = c("HC", "SLE-NP", "SLE-P"), cex = 1.5, fill = color.mixo(1:4), bty = "n")
##


##
# REDUCTION DATA ALL GROPS
list.keepX <- c(5:10,  seq(20, 100, 10))

tune.splsda.srbct <- tune.splsda(SLE_red[,c(2:19)], SLE_red$Phenotype, ncomp = 2, # число компонент берется как число групп-1
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 10)   
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp 
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp] 
#Оптимизированный sPLS-DA
MyResult.splsda <- splsda(SLE_red[,c(2:19)], SLE_red$Phenotype, 
                          ncomp = 2, keepX = select.keepX)
#
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(MyResult.splsda, comp = 1, 
             size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", 
             contrib = "max", legend = FALSE, col.ties="black", ndisplay = 5)
plotLoadings(MyResult.splsda, comp = 2, 
             size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", 
             contrib = "max",ndisplay = 5,  legend = FALSE, col.ties="black")
plotIndiv(MyResult.splsda, ind.names = F, 
          ellipse = T, style = "graphics", abline = TRUE, cex = 2, 
          pch = 19, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, 
          title = "sPLS-DA ordination of samples", size.title = 1.5)
legend("bottomright", legend = c("HC", "SLE-NP", "SLE-P"), cex = 1.5, fill = color.mixo(1:4), bty = "n")



### _____________________________________________________________________________________________

# ALL FEACHURES SLE-NP vs SLE-P
list.keepX <- c(5:10,  seq(20, 100, 10))

tune.splsda.srbct <- tune.splsda(SLE[40:119,c(3:229)], SLE$Phenotype[40:119], ncomp = 2, # число компонент берется как число групп-1
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 10)   
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp 
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp] 
#Оптимизированный sPLS-DA
MyResult.splsda <- splsda(SLE[40:119,c(3:229)], SLE$Phenotype[40:119], 
                          ncomp = 2, keepX = select.keepX)
#
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(MyResult.splsda, comp = 1, 
             size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", 
             contrib = "max", legend = FALSE, col.ties="black", ndisplay = 5)
plotLoadings(MyResult.splsda, comp = 2, 
             size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", 
             contrib = "max",ndisplay = 5,  legend = FALSE, col.ties="black")
plotIndiv(MyResult.splsda, ind.names = F, 
          ellipse = T, style = "graphics", abline = TRUE, cex = 2, 
          pch = 19, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, 
          title = "sPLS-DA ordination of samples", size.title = 1.5)
legend("bottomright", legend = c("SLE-NP", "SLE-P"), cex = 1.5, fill = color.mixo(1:4), bty = "n")
##


##
# REDUCTION DATA SLE-NP vs SLE-P
list.keepX <- c(5:10,  seq(20, 100, 10))

tune.splsda.srbct <- tune.splsda(SLE_red[40:119,c(2:19)], SLE_red$Phenotype[40:119], ncomp = 2, # число компонент берется как число групп-1
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 10)   
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp 
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp] 
#Оптимизированный sPLS-DA
MyResult.splsda <- splsda(SLE_red[40:119,c(2:19)], SLE_red$Phenotype[40:119], 
                          ncomp = 2, keepX = select.keepX)
#
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(MyResult.splsda, comp = 1, 
             size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", 
             contrib = "max", legend = FALSE, col.ties="black", ndisplay = 5)
plotLoadings(MyResult.splsda, comp = 2, 
             size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", 
             contrib = "max",ndisplay = 5,  legend = FALSE, col.ties="black")
plotIndiv(MyResult.splsda, ind.names = F, 
          ellipse = T, style = "graphics", abline = TRUE, cex = 2, 
          pch = 19, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, 
          title = "sPLS-DA ordination of samples", size.title = 1.5)
legend("bottomright", legend = c("SLE-NP", "SLE-P"), cex = 1.5, fill = color.mixo(1:4), bty = "n")



















