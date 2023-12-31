#  Author: Xin Wang  flyboyleo@gmail.com
setwd("~/Desktop/FOXA1_motif_freq/120/") 

library(randomForest)
#dd<-read.table("tmp4.out")[-1]
dd2<-read.table("~/Desktop/FOXA1_motif_freq/120/tmp4.out")[-1,-1]
motifName <- read.table("~/Desktop/FOXA1_motif_freq/120/tmp4.out")[-1,1]

motif<-sapply(as.character(motifName),function(x){
  strsplit(x," ")[[1]][1]
})
# -- --------------------------------------------------------------------

labels<-as.factor(c(rep(1,500),rep(2,500),rep(3,500),rep(4,500)))

dd<-t(dd2)



# tree500 importance-----------------------------------------------------------------

set.seed(123)
randomForest(dd,ntree = 500, y = as.factor(labels),
             importance = T, do.trace = 100, proximity = T,keep.inbag = TRUE) -> RFmodel_500
save(RFmodel_500, file = "RFmodel_120_tree500_time1.Rdata")

#cbind(as.character(motif), RFmodel2_500$importance) -> importance_RF 

#write.table(importance_RF, "importance_RF_120_tree500.xls",sep = "\t", quote = F,col.names = NA)

as.data.frame(cbind(as.character(motif), importance(RFmodel_500,scale = T)),stringsAsFactors = F) -> importance_RF_scaleT1 

table(importance_RF_scaleT1$MeanDecreaseAccuracy>0)


# Feature selection by Mean Decrease Accuracy -----------------------------


dd[,importance_RF_scaleT1$MeanDecreaseAccuracy>0]->dd_filtered

set.seed(321)
randomForest(dd_filtered,ntree = 500, y = as.factor(labels),
             importance = T, do.trace = 100, proximity = T,keep.inbag = TRUE) -> RFmodel2_500
save(RFmodel2_500, file = "RFmodel_120_tree500_time2.Rdata")


# Tuning ------------------------------------------------------------------


a <- tuneRF(dd_filtered, labels)
b <- tuneRF(dd_filtered, labels,improve=0.01,stepFactor = 1.2,mtryStart = 74)
c <- tuneRF(dd_filtered, labels,improve=0.01,stepFactor = 1.2,mtryStart = 88)
d <- tuneRF(dd_filtered, labels,improve=0.01,stepFactor = 1.2,mtryStart = 105)

set.seed(231)
randomForest(dd_filtered,ntree = 500, y = as.factor(labels),
             importance = T, do.trace = 100, proximity = T,keep.inbag = TRUE, mtry = 105) -> RFmodel3_500
save(RFmodel3_500, file = "RFmodel_120_tree500_time3_FOXA1.Rdata")

plot(RFmodel_500$err.rate[,1], ylim = c(0.5,0.8), type = "l")
plot(RFmodel2_500$err.rate[,1], ylim = c(0.5,0.8), type = "l")
plot(RFmodel3_500$err.rate[,1], ylim = c(0.5,0.8), type = "l")

plot(sort(importance(RFmodel_500,scale = T)[,5], decreasing = T), type = "l")
plot(sort(importance(RFmodel2_500,scale = T)[,5], decreasing = T), type = "l")

pdf(file = "importance_curve_FOXA1.pdf", 4, 4)

plot(sort(importance(RFmodel3_500,scale = T)[,5], decreasing = T), type = "l", xlab = "Filtered motifs ranked by MDA"
     , ylab = "Mean Decrease Accuracy value")

dev.off()

importance(RFmodel3_500,scale = T) -> importance_RF_scaleT3
row.names(importance_RF_scaleT3)<-motif[importance_RF_scaleT1$MeanDecreaseAccuracy>0]
                    
importance_RF_scaleT3[importance_RF_scaleT3[,5]>6,1:5]->FOXimportanceBigerthan6 

FOXimportanceBigerthan6[order(FOXimportanceBigerthan6[,5],decreasing = T),]->topMDA

colnames(topMDA) <- c("HEK293T", "HepG2", "K562", "MCF-7", "MDA")

library(gplots)
color <- colorRampPalette(c("blue", "light yellow", "red"))
maintxt<-"Top motifs by MDA of the FOXA1 dataset"
pdf("topMDA_FOXA1.pdf", height=8,width=8)

heatmap.2(topMDA[,1:4],scale="none",Rowv=F,dendrogram="none",Colv=F, #labRow = "", labCol = "", 
          trace="none",sepcolor="black",
          col=color, main = maintxt)

dev.off()


# inTrees ------------------------------------------------------------


library(inTrees)
library(randomForest) 

X <- dd_filtered  # X: predictors
target <- as.factor(labels)  # target: class
#rf <- randomForest(X, as.factor(target))
load("RFmodel_120_tree500_time3_FOXA1.Rdata")
RFmodel3_500->rf
rm(RFmodel3_500)
treeList <- RF2List(rf)  # transform rf object to an inTrees' format
exec <- extractRules(treeList, X)  # R-executable conditions
exec[1:2,]


ruleMetric <- getRuleMetric(exec,X,target)  # get rule metrics
ruleMetric[1:2,]

ruleMetric <- pruneRule(ruleMetric, X, target)
ruleMetric[1:2,]

colnames(X) <- motif[importance_RF_scaleT1$MeanDecreaseAccuracy>0]

ruleMetric <- selectRuleRRF(ruleMetric, X, target)
learner <- buildLearner(ruleMetric, X, target)

motif2 <- motif[importance_RF_scaleT1$MeanDecreaseAccuracy>0]

learner2 <- presentRules(learner, motif2)


readableRules <- presentRules(ruleMetric, motif2)
readableRules[1:2, ]


write.table(learner2, file = "learner2_120_tree500.xls", sep = "\t", col.names = NA)

write.table(readableRules, file = "readableRules_120_tree500.xls", sep = "\t", col.names = NA)





