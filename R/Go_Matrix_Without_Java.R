################################
## Plot Function 
################################
doGoCorrPlot<- function(data,pref = "out"){
        
        
        # data inform of data.frame. input will be output of java program GoMatrix. Manually you have to add last column as label in the input file .
        library(corrplot)
        ## prepare lebel from given input file 
        data <- mat
        ## Do clustering. rearrange matrix row and column before clustering.
        M = cor(data, method="pearson") # do correlation 
        go_clst<-hclust(dist(M)) #do clustering of correlated data
        x<-M[go_clst$order,] # rearrange matrix with clustering order, rearrange row first
        x<-x[,go_clst$order] # rearrange matrix with clustering order, rearrange column now 
        head(colnames(x)) # chk arranged columns 
        head(rownames(x)) # chk arranged rows 
        attributes(go_clst)  
        
        ## generate plot 
        png(file=paste(pref,"MyGoClusterMap.png",sep="_"),height = 5000,width=5000) # create file
        corrplot(as.matrix(x), type="full", method="color",is.corr = F,diag=T,order ="original",tl.cex=1,cl.cex = 5) #do corr plot without  clustering & for single sample 
        #legend("topleft", legend = names(s), col=s, pch="---",cex=5,text.col = s)  # put legend on plot
        dev.off()
        ## write go in order of clustering 
        write.table(rownames(x),file = paste(pref,"MyGoClusterMap.txt",sep="_"),sep = "\t",quote = F)
        
}

################################
## Process input data 
################################

matrixFile="./sampleData/goMatrix.txt" ## input file contains two column. Go and associated genes. Genes MUST BE sepereted by ":"
data <- read.table(matrixFile,sep="\t",quote=NULL, comment='',check.names=F)
head(data)

dat <- sapply(as.character(data[,2]), function(X){return(strsplit(X,":"))}) ## split by deleminator 
names(dat) <- as.character(data[,1]) ## assign names to each elem of list 

mat <- matrix(nrow = length(dat),ncol = length(dat))
colnames(mat) <- names(dat)
row.names(mat) <- names(dat)

combinations <- combn(c(1:length(dat)),2)
numberOfGenesOverlapped  <- apply(combinations, 2, function(X) { 
        row = X[1] 
        col = X[2] 
        return(length(intersect(dat[[row]],dat[[col]])))
        
})

## assign value to matrix.
## chk 
for(i in c(1: ncol(combinations))){
        indx <- combinations[,i]
        mat[indx[1],indx[2]] <- numberOfGenesOverlapped[i]
        mat[indx[2],indx[1]] <- numberOfGenesOverlapped[i]
}

## assign value to diagonal elements 
for(i in c(1: dim(mat)[1])){
        
        mat[i,i] <- 1
}

prefix ="40_updated" ## prefix for output plot 

################################
## Run Function
################################

doGoCorrPlot(mat , pref = prefix)

