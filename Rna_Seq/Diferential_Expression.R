#Rscript

#Analysing Diferential Expresssion
library(edgeR)
setwd("/home/scratch60/gabi_miRNA/align")
d=read.table("allmatrix", header=T, stringsAsFactors=F, row.names=1)
aux=d[, colSums(d) !=0]
dt=t(aux)
input=dt[ , order(colnames(dt))]
grupo=factor(c("controle","controle","controle","controle","controle","controle","controle","controle","controle","controle","controle","controle","pep","pep","pep","pep","pep","pep","pep","pep","pep","pep","pep","pep","pep_tratado","pep_tratado","pep_tratado","pep_tratado","pep_tratado","pep_tratado","pep_tratado","pep_tratado","pep_tratado","pep_tratado","pep_tratado","pep_tratado"))
y = DGEList(counts=input,group=grupo)
keep = rowSums(cpm(y)>1) >= 2
y=y[keep, , keep.lib.sizes=FALSE]
y=calcNormFactors(y)
design <- model.matrix(~grupo)
y=estimateDisp(y,design)
fit=glmQLFit(y,design)
#Comparing 2 vs 1
qlf.2vs1 <- glmQLFTest(fit, coef=2)
#Comparing 3 vs 1
qlf.3vs1 <- glmQLFTest(fit, coef=3)
#Comparing 3 vs 2
qlf.3vs2 <- glmQLFTest(fit, contrast=c(0,-1,1))
out3vs2=topTags(qlf.3vs2, n=nrow(qlf.3vs2$table))
out2vs1=topTags(qlf.2vs1, n=nrow(qlf.2vs1$table))
out3vs1=topTags(qlf.3vs1, n=nrow(qlf.3vs1$table))
table.3vs2=out3vs2$table
table.2vs1=out2vs1$table
table.3vs1=out3vs1$table
filt.st.3vs2=table.3vs2[table.3vs2$FDR < 0.01 & abs(table.3vs2$logFC) > 2 ,]
filt.nst.3vs2=table.3vs2[table.3vs2$FDR < 0.05 & abs(table.3vs2$logFC) > 1.5 ,]
filt.st.3vs1=table.3vs1[table.3vs1$FDR < 0.01 & abs(table.3vs1$logFC) > 2 ,]
filt.nst.3vs1=table.3vs1[table.3vs1$FDR < 0.05 & abs(table.3vs1$logFC) > 1.5 ,]
filt.st.2vs1=table.2vs1[table.2vs1$FDR < 0.01 & abs(table.2vs1$logFC) > 2 ,]
filt.nst.2vs1=table.2vs1[table.2vs1$FDR < 0.05 & abs(table.2vs1$logFC) > 1.5 ,]

#Building Volcano plots
ggplot(table.2vs1, aes(x=table.2vs1$logFC, y=-log2(table.2vs1$FDR))) + geom_point(size = 2.5, alpha = 1, na.rm = T) + xlab("logFC") + ylab("-log(FDR)") + geom_hline(yintercept = 1.3, linetype="dashed") + geom_vline(xintercept=3.8, linetype="dashed") + geom_vline(xintercept=-3.8, linetype="dashed")


#heatmap for differentially expressed genes

teste1=rownames(filt.nst.3vs2)
teste2=rownames(filt.nst.2vs1)
teste3=rownames(filt.nst.3vs1)
diffrows=unique(c(teste1,teste2,teste3))
heatmap(dt[diffrows,])

#Calculating number of miRNAs expressed per sample

for (i in colnames(dt)){
  col=dt[,i]
  d=col[col>=1]
  print(i)
  print (length(d)) }

#Calculating number of miRNA reads per sample

apply(dt,2,sum)

