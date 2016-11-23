################################
## load  Function 
################################
doGoCorrPlot<- function(filePath = "./sampleData/goMatrix.txt", pref = "out", outImageHeight = 4000, outImageWidth= 4000,colorLabelSize = 3 ,GOLabelSize =1){
        
        
        # install required packages 
        
        if(!require(corrplot)){
                install.packages("corrplot")
        }
        library(corrplot)
        
        ################################
        ## Process input data 
        ################################
      
        matrixFile=filePath ## input file contains two column. Go and associated genes. Genes MUST BE sepereted by ":"
        data <- read.table(matrixFile,sep="\t",quote=NULL, comment='',check.names=F)
        head(data)
        
        dat <- sapply(as.character(data[,2]), function(X){return(strsplit(X,":"))}) ## split by deleminator 
        names(dat) <- as.character(data[,1]) ## assign names to each elem of list 
        
        mat <- matrix(nrow = length(dat),ncol = length(dat)) # create empty matrix
        colnames(mat) <- names(dat) # assign colnames 
        row.names(mat) <- names(dat) # assign row names
        
        combinations <- combn(c(1:length(dat)),2) ### create all possible combination for given dimension of matrix 
        
        # get number of overlapped genes between each pair of GO 
        numberOfGenesOverlapped  <- apply(combinations, 2, function(X) { 
                row = X[1] 
                col = X[2] 
                return(length(intersect(dat[[row]],dat[[col]])))
                
        })
        
        ## assign value to matrix.
        
        for(i in c(1: ncol(combinations))){
                indx <- combinations[,i]
                mat[indx[1],indx[2]] <- numberOfGenesOverlapped[i]
                mat[indx[2],indx[1]] <- numberOfGenesOverlapped[i]
        }
        
        ## assign value to diagonal elements 
        for(i in c(1: dim(mat)[1])){
                
                mat[i,i] <- 1
        }
        
        ## Do clustering. rearrange matrix row and column before clustering.
        M = cor(mat, method="pearson") # do correlation 
        go_clst<-hclust(dist(M)) #do clustering of correlated data
        x<-M[go_clst$order,] # rearrange matrix with clustering order, rearrange row first
        x<-x[,go_clst$order] # rearrange matrix with clustering order, rearrange column now 
        head(colnames(x)) # chk arranged columns 
        head(rownames(x)) # chk arranged rows 
        attributes(go_clst)  
        
        ## generate plot 
        png(file=paste(pref,"MyGoClusterMap.png",sep="_"),height = outImageHeight,width=outImageHeight) # create file
        corrplot(as.matrix(x), type="full", method="color",is.corr = F,diag=T,order ="original",tl.cex=GOLabelSize,cl.cex = colorLabelSize) #do corr plot without  clustering & for single sample 
        #legend("topleft", legend = names(s), col=s, pch="---",cex=5,text.col = s)  # put legend on plot
        dev.off()
        ## write go in order of clustering 
        write.table(data[go_clst$order,],file = paste(pref,"MyGoClusterMap.txt",sep="_"),sep = "\t",quote = F,row.names = F,col.names = c("GO_In_Order_Of_Clustering","Genes"))
        
}

################################
## Run Function
################################
prefix ="Out_GO_Matrix" ## prefix for output plot 

doGoCorrPlot(filePath = "./sampleData/goMatrix.txt", pref = prefix)

