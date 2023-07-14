load("~/Desktop/FOXA1_motif_freq/120/RFmodel_120_tree1k_FOXA1.Rdata")
dd1<-read.table("~/Desktop/mystery/mystery1.bed.500.txt.txt.motifs")[-1,-1]
dd2<-read.table("~/Desktop/mystery/mystery2.bed.500.txt.txt.motifs")[-1,-1]
dd3<-read.table("~/Desktop/mystery/mystery3.bed.500.txt.txt.motifs")[-1,-1]
data1<-t(dd1)
data2<-t(dd2)
data3<-t(dd3)
predict(RFmodel_1000, data1, type="response")
predict(RFmodel_1000, data2, type="response")
predict(RFmodel_1000, data3, type="response")

# I am not confident about those results, since the regions
# from the same cell type seem not to be assigned to the 
# same class according to this model. 
# This might imply that the model is not accurate enough, e.g.,
# due to overfitting, or those ccREs have different motifs
# responsible for FOXA1 binding. 


