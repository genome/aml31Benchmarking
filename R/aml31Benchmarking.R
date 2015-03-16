library(venneuler)


addKey <- function(df){
    df$key = paste(df[,1],df[,2],df[,4],df[,5],sep="_")
    return(df)
}

##----------------------------------------------
subsetByVaf <- function(df,range){
  if("tum_vaf" %in% names(df)){
    return(df[df$tum_vaf >= range[1] & df$tum_vaf <= range[2],])
  }
  return(df)
}

##----------------------------------------------
##take two dataframes, each with a "key" column, and create a proportional venn diagram of their overlap
createVenn <- function(a, b, range, label){

  #use range later to restarict to certain vafs
  a = subsetByVaf(a,range)
  b = subsetByVaf(b,range)

  allkeys = unique(c(as.character(a$key), as.character(b$key)))
  keyFound = cbind(allkeys %in% a$key, allkeys %in% b$key)
  A = length(which(keyFound[,1]==1 & keyFound[,2]==0))
  B = length(which(keyFound[,1]==0 & keyFound[,2]==1))
  AB = length(which(keyFound[,1]==1 & keyFound[,2]==1))

  v <- venneuler(c(A=A, B=B, "A&B"=AB))
  plot(v,main=paste("Comparison to",label,"list"))

  col.fn <- function(col, alpha=0.3) {
    col <- hcl(col * 360, 130, 60)
    col <- col2rgb(col)/255
    col <- rgb(col[1, ], col[2, ], col[3, ], alpha)
    col
  }
  vcol <- c(col.fn(v$colors),rgb(130/255,209/255,201/255))
  vlabs <- v$labels

  legend(0.05, .9, legend=c(paste("Truth-set Only:",A), paste("Uploaded Only:",B), paste("Shared:",AB)), fill=vcol, x="topleft")
  #print sens/spec in bottom left
  usr <- par( "usr" )
  text( usr[ 1 ], usr[ 3 ]+0.1, paste("Sensitivity:",round((AB/A),4)), adj = c( 1, 0 ), pos=4)
  text( usr[ 1 ], usr[ 3 ]+0.05, paste("Positive Predictive Value:",round((AB/(AB+B)),4)), adj = c( 1, 0 ), pos=4)
  text( usr[ 1 ], usr[ 3 ]+0.15, paste("VAF range: ",range[1],"-",range[2],sep=""),adj = c( 1, 0 ), pos=4)
}

##----------------------------------------------
#create histogram showing performance per VAF
createHist <- function(b, a, label){
  v = NULL;  
  for(vaf in seq(0,95,5)){
    av = a[a$tum_vaf > vaf & a$tum_vaf <= vaf+5,]
    allkeys = unique(c(as.character(av$key), as.character(b$key)))
    keyFound = cbind(allkeys %in% av$key, allkeys %in% b$key)
    A = length(which(keyFound[,1]==1 & keyFound[,2]==0))
    AB = length(which(keyFound[,1]==1 & keyFound[,2]==1))

    v = rbind(v, data.frame(vaf=vaf, val=A, shared=AB))
  }
  print(v)
  #cat(v$val,file=stderr())
  vals = rep(v$vaf,v$val)
  hist(vals,breaks=seq(0,95,5),col=rgb(0,1,0,0.3),xlab="Tumor VAF",main="Valid uploaded variants by VAF")
  mtext(paste("comparison to",label,"list"),cex=0.9)
  vals = rep(v$vaf,v$shared)
  hist(vals,breaks=seq(0,95,5),col=rgb(0,0,1,0.3),add=T,xlab="",main="")
  legend("topright",fill=c(rgb(0,1,0,0.3),rgb(0,0,1,0.3)), legend=c("Truth-set Variants","Uploaded Variants in Truth-set"))
}

##----------------------------------------------
#strip the 'key' column from a data frame for display purposes
stripKeyCol <- function(df){
  return(df[,!(names(df) %in% c("key"))])
}

##----------------------------------------------
#label uploaded variants as found in truthset or not
labelFound<-function(upload,truth){
  upload$in_validation = as.character(upload$key) %in% as.character(truth$key)
  return(stripKeyCol(upload))
}


##---------------------------------------------------------
#truthList

benchmark <- function(variantFile,tablePrefix,plotFile){

  data(list.platinum)
  data(list.gold)

  pdf(plotFile,width=8,height=5)

  #load("../data/list.platinum.Rdata")
  userData = read.table(variantFile,header=F,sep="\t", col.names=c("Chr","St","Sp","Ref","Var"))
  createVenn(addKey(userData), addKey(list.platinum), range=c(0,50), label="Platinum")
  createHist(addKey(userData), addKey(list.platinum), label="Platinum")
  userData = labelFound(addKey(userData), addKey(list.platinum))
  write.table(userData,paste(tablePrefix,".platinum",sep=""),sep="\t",row.names=F,quote=F)

  #load("../data/list.gold.Rdata")
  userData = read.table(variantFile,header=F,sep="\t", col.names=c("Chr","St","Sp","Ref","Var"))
  createVenn(addKey(userData), addKey(list.gold), range=c(0,50), label="Gold")
  createHist(addKey(userData), addKey(list.gold), label="Gold")
  userData = labelFound(addKey(userData), addKey(list.gold))
  write.table(userData,paste(tablePrefix,".gold",sep=""),sep="\t",row.names=F,quote=F)
  
  dev.off()
}
