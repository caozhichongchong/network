################## C-score calculation and correlation analysis ###############################
library(vegan)
library(igraph)
library(Hmisc)

co_occurrence_network_cyto<-function(matrix,alpha,p.cutoff,Direction=1){
  #C-score calculation using null model
  matrix1<-matrix
  matrix1[matrix1>0]<-1
  
  ### creating node file
  node<-matrix(0,ncol=1,nrow=nrow(matrix))
  rownames(node)<-rownames(matrix)
  node[,1]<-rowMeans(matrix)
  
  
  #co_occurrence<-oecosimu(t(matrix1),nestedchecker,"swap",nsimul=99)
  
  #correlation analysis based on spearman's co-efficient
  matrix.dist<-rcorr(t(matrix),type="spearman")
  matrix.cor<-matrix.dist$r
  matrix.cor.p<-matrix.dist$P
  
  #Multiple testing correction using Benjamini-Hochberg standard false discovery rate correction ("FDR-BH")
  matrix.cor.p <- p.adjust(matrix.cor.p, method="BH")
  
  ##1.Consider positive cooccurence at given coefficient (alpha) and p-value cutoffs
  if(Direction==1)
  {matrix.cor1<-matrix.cor
  matrix.cor1.p<-matrix.cor.p
  matrix.cor1[which(matrix.cor1 <= alpha)]=0
  matrix.cor1[which(matrix.cor1.p>p.cutoff)]=0
  matrix.cor1[which(is.na(matrix.cor1))]=0
  ## delete those rows and columns with sum = 0
  matrix.cor1<-matrix.cor1[which(rowSums(matrix.cor1)!=0),]
  matrix.cor1<-matrix.cor1[,which(colSums(matrix.cor1)!=0)]}
  
  ##2.Consider netagive cooccurence at given coefficient (-alpha) and p-value cutoffs
  else if(Direction==2)
  {matrix.cor2<-matrix.cor
  matrix.cor2.p<-matrix.cor.p
  matrix.cor2[which(matrix.cor2 > (-alpha))]=0
  matrix.cor2[which(matrix.cor2.p>p.cutoff)]=0
  matrix.cor3[which(is.na(matrix.cor2))]=0
  ## delete those rows and columns with sum = 0
  matrix.cor2<-matrix.cor2[which(rowSums(matrix.cor2)!=0),]
  matrix.cor2<-matrix.cor2[,which(colSums(matrix.cor2)!=0)]}
  
  ##3.Consider both positive and netagive cooccurence at given coefficient (alpha) and p-value cutoffs
   else if(Direction==3)
  {matrix.cor3<-matrix.cor
   matrix.cor3.p<-matrix.cor.p
   matrix.cor3[which(matrix.cor3>=(-alpha) & matrix.cor3 <= alpha)]=0
   matrix.cor3[which(matrix.cor3.p>p.cutoff)]=0
   matrix.cor3[which(is.na(matrix.cor3))]=0
  ## delete those rows and columns with sum = 0
   matrix.cor3<-matrix.cor3[which(rowSums(matrix.cor3)!=0),]
   matrix.cor3<-matrix.cor3[,which(colSums(matrix.cor3)!=0)]}
  
###create edge file
if(Direction==1)
{edge<-matrix(0,ncol=3,nrow=(nrow(matrix.cor1)*(nrow(matrix.cor1)-1)))
k<-1
for(i in 1:(nrow(matrix.cor1)))
{
  for(j in i:nrow(matrix.cor1))
  {
    if(matrix.cor1[i,j]!=0 && matrix.cor1[i,j]!=1)
    {
      edge[k,1]<-rownames(matrix.cor1)[i]
      edge[k,2]<-rownames(matrix.cor1)[j]
      edge[k,3]<-as.numeric(matrix.cor1[i,j])
      k<-k+1
    }
  }
}}
  
else if(Direction==2)
  {edge<-matrix(0,ncol=3,nrow=(nrow(matrix.cor2)*(nrow(matrix.cor2)-1)))
  k<-1
  for(i in 1:(nrow(matrix.cor2)))
  {
    for(j in i:nrow(matrix.cor2))
    {
      if(matrix.cor2[i,j]!=0 && matrix.cor2[i,j]!=1)
      {
        edge[k,1]<-rownames(matrix.cor2)[i]
        edge[k,2]<-rownames(matrix.cor2)[j]
        edge[k,3]<-as.numeric(matrix.cor2[i,j])
        k<-k+1
      }
    }
  }}
  
else if(Direction==3)
  {edge<-matrix(0,ncol=3,nrow=(nrow(matrix.cor3)*(nrow(matrix.cor3)-1)))
  k<-1
  for(i in 1:(nrow(matrix.cor3)))
  {
    for(j in i:nrow(matrix.cor3))
    {
      if(matrix.cor3[i,j]!=0 && matrix.cor3[i,j]!=1)
      {
        edge[k,1]<-rownames(matrix.cor3)[i]
        edge[k,2]<-rownames(matrix.cor3)[j]
        edge[k,3]<-as.numeric(matrix.cor3[i,j])
        k<-k+1
      }
    }
  }}
  # append the output into results
 result<-list()
  #result$co_occurrence<-co_occurrence
  result$node<-node
  result$edge<-edge[1:(k-1),]
  
# result$matrix.cor1<-matrix.cor1
#  result$graph1<-g1
  
#  result$matrix.cor2<-matrix.cor2
 # result$graph2<-g2
  
#result$matrix.cor3<-matrix.cor3
#  result$graph3<-g3
  return(result)
} 
