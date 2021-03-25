library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
options(bitmapType='cairo')
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cbPalette <- c("#999999","#F0E442","#E69F00","#D55E00")

parameters=read.table("R_parameters.txt",sep="\t",header=FALSE,row.names=1,stringsAsFactors = FALSE)
sample_list=read.table("limma_sample_list.txt",sep="\t",header=TRUE,row.names=1)
data=read.table("count_file.txt",sep="\t",header=TRUE,row.names=1)

capabilities()

print(parameters["contrasts",1])




group=factor(sample_list[,1])

design <- model.matrix(~ 0+group)
colnames(design)=levels(group)
print(design)
filterCount=strtoi(parameters["filter_count",1])
filterCountType=parameters["filter_count_type",1]

if (filterCountType=="raw_counts"){
  filteredData=subset(data,rowSums(data)>=filterCount)
  dge=DGEList(counts=filteredData,group=group)

} else if (filterCountType=="CPM"){
    dge=DGEList(counts=data,group=group)
    cutoff<-as.vector(cpm(filterCount,mean(dge$samples$lib.size) ) )
    keep<-rowSums( cpm(dge) >cutoff ) >=min(table(sample_list[,1]) )
    dge<-dge[keep,]
}

geneCount=nrow(dge$counts)
print(paste("GeneCount:",geneCount,sep=""))

command=paste("makeContrasts(",parameters["contrasts",1],",levels=design)",sep="")
print("command")
print(command)
contrast.matrix = eval(parse(text=command))

dge <- calcNormFactors(dge,method=parameters["calc_norm_factors_method",1])
par(mar=c(1,1,1,1))
png(filename="mean_variance.png",height=1500,width=1500,res=300)
v <- voom(dge,design,plot=TRUE)
dev.off()
fit <- lmFit(v,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)



if (nrow(sample_list)>=3){
  mds=plotMDS(v)
  dev.off()
  mdsDF=data.frame(mds$x,mds$y,rownames(sample_list))
   	
	rownames(mdsDF)=rownames(sample_list)
	colnames(mdsDF)=c("x","y","Name")
  png(filename="mds_plot.png",height=1500,width=1500,res=300)

  gp=ggplot(mdsDF, aes(x,y )) + geom_point( size = 1.5, color="gray") +  geom_text_repel(aes(label = Name), point.padding = unit(0.25, "lines"), box.padding = unit(0.25, "lines"), nudge_y = 0.1) + theme_bw(base_size = 10)
	print(gp)

  dev.off()
}

sigfig <- function(vec, n=5){ 
    formatC(signif(vec,digits=n), digits=n,format="fg", flag="#") 
}

contrastList=unlist(strsplit(parameters["contrasts",1],","))

contrast_files=NULL

for (con in contrastList){
  pValueAdjust=parameters["p_value_adjust",1]
  
  nameCon=unlist(strsplit(con,"="))
  contrastName=nameCon[1]
  c=nameCon[2]

  print("coef info")
  print(contrastName)
  print(c)

  resultTable=topTable(fit3, coef=contrastName, adjust=pValueAdjust,sort="none",n=geneCount)
  topHits=topTable(fit3, coef=contrastName, adjust=pValueAdjust,sort="P",n=25)
  sortedResultTable=data.frame(topTable(fit3, coef=contrastName, adjust=pValueAdjust,sort="P",n=geneCount))
  
  sortedResultTable$Name=""
  sortedResultTable$Gene=rownames(sortedResultTable)
  sortedResultTable[1:min(c(25,nrow(sortedResultTable))),]$Name=as.character(sortedResultTable[1:min(c(25,nrow(sortedResultTable))),]$Gene)
  sortedResultTable$AdjustedPvalue=">=0.1"
  print("test .1")
  if (nrow(sortedResultTable[sortedResultTable$adj.P.Val<.1,])>0){
    print (nrow(sortedResultTable[sortedResultTable$adj.P.Val<.1,]))
    sortedResultTable[sortedResultTable$adj.P.Val<.1,]$AdjustedPvalue="<0.1"
  }
  print("test .05")
  if (nrow(sortedResultTable[sortedResultTable$adj.P.Val<.05,])>0){
    sortedResultTable[sortedResultTable$adj.P.Val<.05,]$AdjustedPvalue="<0.05"
  }
  print("test .01")
  if (nrow(sortedResultTable[sortedResultTable$adj.P.Val<.01,])>0){
    sortedResultTable[sortedResultTable$adj.P.Val<.01,]$AdjustedPvalue="<0.01"
  }
  print("end test")
  #sortedResultTable$AdjustedPvalue=factor(sortedResultTable$AdjustedPvalue)
  
  vals=(c(">=0.1","<0.1","<0.05","<0.01"))
  #vals=(c("<0.01","<0.05","<0.1",">=0.1"))
  sortedResultTable$AdjustedPvalue=factor(sortedResultTable$AdjustedPvalue,levels=vals)
  
  if (pValueAdjust=='fdr'){
      pValueAdjust="FDR"
  }
  
  colnames(sortedResultTable)[9]=pValueAdjust
  print("MA PLot")
  print(colnames(sortedResultTable))
  graphics.off()
  
  png(filename=paste("MA_plot",contrastName,"png",sep="."),height=1500,width=2000,res=300)
  gp=ggplot(sortedResultTable, aes(AveExpr, logFC)) + geom_point(aes_string(fill = pValueAdjust),size = 1.5, color="black",pch=21) +  geom_text_repel(aes(label = Name), point.padding = unit(0.25, "lines"), box.padding = unit(0.25, "lines"), nudge_y = 0.1) + theme_bw(base_size = 10)+scale_fill_manual(values=cbPalette)
  print(gp)

  dev.off()
  
  print("Volcano")
  png(filename=paste("volcano_plot",contrastName,"png",sep="."),height=1500,width=1500,res=300)
  sortedResultTable$mLog10=-log10(as.numeric(sortedResultTable$adj.P.Val))
  gp=ggplot(sortedResultTable, aes(logFC,mLog10 )) + geom_point( size = 1.5, color="gray") +  geom_text_repel(aes(label = Name), point.padding = unit(0.25, "lines"), box.padding = unit(0.25, "lines"), nudge_y = 0.1) + theme_bw(base_size = 10)+labs(y = paste("-log10(",pValueAdjust,")",sep=""))+ylim(c(0,max(sortedResultTable$mLog10)*1.2))
  print(gp)
  dev.off()
  
  sortedResultTable$PLog10=-log10(as.numeric(sortedResultTable$P))
  png(filename=paste("volcano_plot.pvalue",contrastName,"png",sep="."),height=1500,width=1500,res=300)
 
  gp=ggplot(sortedResultTable, aes(logFC,PLog10 )) + geom_point( size = 1.5, color="gray") +  geom_text_repel(aes(label = Name), point.padding = unit(0.25, "lines"), box.padding = unit(0.25, "lines"), nudge_y = 0.1) + theme_bw(base_size = 10)+labs(y = "-log10(pValue)")+ylim(c(0,max(sortedResultTable$PLog10)*1.2))
  print("C")
  print(gp)
  print("D")
  dev.off()
  
  print("TABLES")

  vE=data.frame(t(apply(v$E,1,function(x) sigfig(x,5))))
  resultsF=data.frame(t(apply(resultTable,1,function(x) sigfig(x,5))))
  resultsF$P.Value=as.numeric(as.character(resultsF$P.Value))
  resultsF$adj.P.Val=as.numeric(as.character(resultsF$adj.P.Val))

  cve=colnames(vE)

  colnames(vE)=paste(cve,"_log2CPM",sep="")

  finalResultTable=cbind(rownames(resultTable),resultsF,vE)

  colnames(finalResultTable)[1]="gene"
  outTable <- finalResultTable[order(finalResultTable$P.Value),] 
  contrastFile=paste("results",contrastName,"txt",sep=".")
  write.table(outTable,file=contrastFile,sep="\t",quote=F,row.names=F)
  contrast_files=rbind(contrast_files,data.frame(contrastName))
    
  
  gseaTable=cbind(rownames(vE),rownames(vE),vE)
  colnames(gseaTable)[1]="Name"
  colnames(gseaTable)[2]="Description"
  write.table(gseaTable,file=paste("GSEA.input",contrastName,"txt",sep="."),sep="\t",quote=F,row.names=F)
  
  tStat=data.frame(gseaTable$Name,resultsF$t)
  write.table(tStat,file=paste("GSEA.tStat",contrastName,"txt",sep="."),sep="\t",quote=F,row.names=F)
}
write.table(contrast_files,file="contrast_files.txt",sep="\t",quote=F,row.names=F,col.names = F)
