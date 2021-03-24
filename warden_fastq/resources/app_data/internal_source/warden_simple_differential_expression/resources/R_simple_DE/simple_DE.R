library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
options(bitmapType='cairo')
writeLines(capture.output(sessionInfo()), "sessionInfo_simple.txt")
cbPalette <- c("#999999","#F0E442","#E69F00","#D55E00")

print("reading parameters")
parameters=read.table("R_parameters_simple.txt",sep="\t",header=FALSE,row.names=1,stringsAsFactors = FALSE)
print("reading sample_list")
sample_list=read.table("sample_list.txt",sep="\t",header=TRUE,row.names=1)
print("Reading count file")
data=read.table("count_file.txt",sep="\t",header=TRUE,row.names=1)


group=factor(sample_list[,1])
design <- model.matrix(~ 0+group)
colnames(design)=levels(group)
print(design)

cs=colSums(data)
norm=t(apply(data,1,function(x) log2(((x+.5)/cs*(10**6)))))
print("MDS")
if (nrow(sample_list)>=3){
  rm=rowMeans(norm)
  ssN=subset(norm,rm>-2)
  mds <- plotMDS(ssN)
  dev.off()
  mdsDF=data.frame(mds$x,mds$y,rownames(sample_list))
   	
	rownames(mdsDF)=rownames(sample_list)
	colnames(mdsDF)=c("x","y","Name")
	print("XX")
  png(filename="mds_plot.norm_CPM.png",height=1500,width=1500,res=300)

  gp=ggplot(mdsDF, aes(x,y )) + geom_point( size = 1.5, color="gray") +  geom_text_repel(aes(label = Name), point.padding = unit(0.25, "lines"), box.padding = unit(0.25, "lines"), nudge_y = 0.1) + theme_bw(base_size = 10)
	print(gp)

  dev.off()
  print("YY")
}

 print("WRITE CPM") 
 sigfig <- function(vec, n=5){ 
    formatC(signif(vec,digits=n), digits=n,format="fg", flag="#") 
 }  
normPretty=data.frame(t(apply(norm,1,function(x) sigfig(x,5))))


write.table(normPretty,file="CPM_normalized_values.txt",sep="\t",quote=F,row.names=T,col.names = NA)


contrast_files=NULL


contrasts=parameters["contrasts",1]

print(contrasts)

if (!is.na(contrasts)) { 

contrastList=unlist(strsplit(contrasts,","))

  for (con in contrastList){
    nameCon=unlist(strsplit(con,"="))
    contrastName=nameCon[1]
    c=nameCon[2]
    cx=unlist(strsplit(c,"-"))
    print(cx[1])
    print(cx[2])

    d=data.frame(design)
    
    c1Cols=as.numeric(row.names(d[d[,cx[1]]==1,]))
    c2Cols=as.numeric(row.names(d[d[,cx[2]]==1,]))
    if (length(c1Cols)==1){
      c1=norm[,c1Cols[1]]
    } else{
      c1=rowMeans(norm[,c1Cols])
    }
    if (length(c2Cols)==1){
      c2=norm[,c2Cols[1]]
    } else{
      c2=rowMeans(norm[,c2Cols])
    }

    df=data.frame(c1,c2)
    
    df$Avg=(df[,1]+df[,2])/2
    df$logFC=(df[,1]-df[,2])
    colnames(df)=c(cx[1],cx[2],"AveExpr","logFC")
    dfPretty=data.frame(t(apply(df,1,function(x) sigfig(x,5))))
    geneNames=rownames(dfPretty)
    dfPretty2=cbind(geneNames,dfPretty)
    dfPretty2$P.value=1
    dfPretty2$t.value=1
    colnames(dfPretty2)[1]="gene"
    
    png(filename=paste("MA_plot","simple_DE",contrastName,"png",sep="."),height=1500,width=2000,res=300)
    gp=ggplot(df, aes(AveExpr, logFC)) + geom_point(size = 1.5, color="black",pch=21)   + theme_bw(base_size = 10)+scale_fill_manual(values=cbPalette)
    print(gp)
    
    dev.off()
    
    write.table(dfPretty2,file=paste("simple_DE",contrastName,"txt",sep="."),sep="\t",quote=F,row.names=FALSE)
    
    contrastFile=paste("simple_DE",contrastName,"txt",sep=".")

    contrast_files=rbind(contrast_files,data.frame(contrastName))
  }
  
  
}
write.table(contrast_files,file="contrast_files.txt",sep="\t",quote=F,row.names=F,col.names = F)
