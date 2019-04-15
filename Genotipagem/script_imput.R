## step 1 ###################

library(dplyr)
bim <- read.table("ptsd_QC_ibd_snpid", header=F)
head(bim)

# novo cabeçalho para merge
colnames(bim) <- c("V1", "SNP", "FID", "V2", "A1","A2")
head(bim)

# quebrar o arquivo referencia em arquivos menores #cria xaa ~ xaf
system(paste("split -l 15000000 cutcols.vcf"))

# Read e merge bim e referência, no fim eu terei 5 arquivos mergidos, cada um com cerca de 60K SNPs
xaa <- read.table("xaa", header=T, sep=",")
head(xaa)
colnames(xaa) <- c("SNP", "REF")
xaabim <- inner_join(bim, xaa, by = "SNP", copy = FALSE)
xaabim <- as.data.frame(xaabim)

xab <- read.table("xab", header=F, sep=",")
colnames(xab) <- c("SNP", "REF")
xabbim <- inner_join(bim, xab, by = "SNP", copy = FALSE)
xabbim <- as.data.frame(xabbim)

xac <- read.table("xac", header=F, sep=",")
colnames(xac) <- c("SNP", "REF")
xacbim <- inner_join(bim, xac, by = "SNP", copy = FALSE)
xacbim <- as.data.frame(xacbim)

xad <- read.table("xad", header=F, sep=",")
colnames(xad) <- c("SNP", "REF")
xadbim <- inner_join(bim, xad, by = "SNP", copy = FALSE)
xadbim <- as.data.frame(xadbim)

xae <- read.table("xae", header=T, sep=",")
colnames(xae) <- c("SNP", "REF")
xaebim <- inner_join(bim, xae, by = "SNP", copy = FALSE)
xaebim <- as.data.frame(xaebim)

xaf <- read.table("xaf", header=T, sep=",")
colnames(xaf) <- c("SNP", "REF")
xafbim <- inner_join(bim, xaf, by = "SNP", copy = FALSE)
xafbim <- as.data.frame(xafbim)

xbim_fim <- rbind(xaabim, xabbim, xacbim, xadbim, xaebim, xafbim)
xbim_fim <- as.data.frame(xbim_fim)
head(xbim_fim)

REF <-xbim_fim[!duplicated(xbim_fim[,2]), ]

#só o que eu preciso, SNP name do bim + REF allele
out <- REF[,c("SNP","REF")]
head(out)

write.table(out, file="xbim_refPTSD.txt", sep=" ", col.names=F, row.names=F, quote=F)


######## step 3 ############


#Pegar RS que deram warning no passo da linha de comando
warning <- read.csv("ptsd_ref.log", header=F, sep=" ")
war <- warning[-(1:23),]
war <- as.data.frame(war)
war[which(1:nrow(war) %% 2 == 0) , ] -> warrsteste
warrsnp <- warrsteste[,-1]
warrsnp$V3 = NULL
warrsnp$V4 = NULL
warrsnp$V5 = NULL
warrsnp$V6 = NULL
warrsnp <- warrsnp[-(217068:217071),] #números de linhas específico do dataset, não copiar; excluir 4 ultimas linhas
warrsnp <- as.data.frame(warrsnp)

write.table(warrsnp, file="war.csv", sep=" ", col.names=F, row.names=F, quote=F)

warcsv <- read.csv("war.csv", header=F, sep=' ')
warrs. <- gsub("[.]$", "", warcsv$V1)
warrs. <- as.data.frame(warrs.)
head (warrs.)
write.table(warrs., file="warrsPTSD.txt", col.names=F, row.names=F, quote=F)


