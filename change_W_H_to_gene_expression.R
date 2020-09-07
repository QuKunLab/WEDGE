
dir = "the_path_of_WEDGE_results"

Gene_name <- read.csv(paste0(dir,"geneName.csv"),header = F)
Cell_name <- read.csv(paste0(dir, "cellName.csv"),header = F)

W<-read.csv(paste0(dir,"W.csv"),header = F)
H<-read.csv(paste0(dir,"H.csv"),header = F)

A<-as.matrix(W)%*%as.matrix(H)
#A<-A%*%as(diag(10000/colSums(A)), "dgCMatrix")
A<-as.matrix(A)
row.names(A) = Gene_name$V1
colnames(A) = Cell_name$V1
