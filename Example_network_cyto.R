setwd("yourdir")
# set up cutoffs
Spearman_cutoff=0.8
P_cutoff=0.01
# input data
Abu=read.table('your_abundance_file.txt',header=T)
Abu<-as.matrix(Abu)
# remove OTUs with total abundance of 0
table<-Abu
table[table>0]<-1
table.generalist<-Abu[which(rowSums(table)>=1),]
Abu<-table.generalist
# calculate spearman value
data<-co_occurrence_network_cyto(Abu,Spearman_cutoff,P_cutoff)
# output node and edge table
nodefile<-data$node
edgefile<-data$edge
write.table(nodefile,"nodetable.txt", col.names = T, sep='\t', quote = FALSE)
write.table(edgefile,"edgetable.txt", col.names = T, sep='\t', quote = FALSE)
