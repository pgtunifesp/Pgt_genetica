########################################################################
# Title: Illumina EPIC methylation array preprocessing, normalization and statistics
# Author: Leticia Spindola
# Date: 23/05/2018
# Purpose: Fazer o QC, normalization, and etc with INPD samples and CAG controls
#######################################################################

library(wateRmelon)
library(minfi)
require(dplyr)
library(ggplot2)
library(qqman)
library(DMRcate)
library(limma)
library(RColorBrewer)
library(Gviz)

dataDirectory <- "~/Documents/ILS_INPD_48"
setwd("~/Documents/ILS_INPD_48")

############################################################
################ IMPORTANDO DADOS EPIC  ####################
############################################################
# Importando os dados fenotipicos das amostras
targets_HRC <- read.csv2(file = "~/Documents/ILS_INPD_48/pData_INPDmeth_final_03_07_2018_minfi.csv", header = T)
class(targets_HRC) # [1] "data.frame"

# Importando os arquivos brutos do EPIC array (.dat)
RGset_HRC <- read.metharray.exp(targets = targets_HRC, extended = TRUE)
dim(pData(RGset_HRC)) # [1] 48 49

annotation(RGset_HRC)
# array                                           annotation 
# "IlluminaHumanMethylationEPIC"                 "ilm10b2.hg19"

# Renomeando os dados brutos com os IDs corretos das amostras
sampleNames(RGset_HRC) <- targets_HRC$sampleID
RGset_HRC
# class: RGChannelSetExtended 
# dim: 1,051,943 48 
# rownames(1051943): 1600101 1600111 ... 99810990 99810992
# colnames(48): I211_W0 I331_W0 ... I236_W1 I378_W1
# colData names(49): BARCODES sampleID ... Basename filenames
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b2.hg19

############################################################
############ IMPORTANDO DADOS CAG CONTROLS  ################
############################################################
targets_controls <- read.csv2(file = "~/Documents/ILS_INPD_48/Controls_CAG_minfi_04_07_2018.csv", header = T)
class(targets_controls) # [1] "data.frame"

RGset_controls <- read.metharray.exp(targets = targets_controls, extended = TRUE)
dim(pData(RGset_controls)) # [1] 140  11

annotation(RGset_controls)
# array                                           annotation 
# "IlluminaHumanMethylation450k"                  "ilmn12.hg19" 

# Renomeando os dados brutos com os IDs corretos das amostras
targets_controls$number <- 1:140
sampleNames(RGset_controls) <- targets_controls$number
RGset_controls
# class: RGChannelSetExtended 
# dim: 622399 140 
# rownames(622399): 10600313 10600322 ... 74810490 74810492
# colnames(140): 1 2 ... 139 140
# Annotation
# array: IlluminaHumanMethylation450k
# annotation: ilmn12.hg19

############################################################
################### CONVERTING ARRAY  ######################
############################################################
# Transformando todas as colunas factor em chacteres
factor_vars <- lapply(colData(RGset_HRC), class) == "factor"
colData(RGset_HRC)[, factor_vars] <- lapply(colData(RGset_HRC)[, factor_vars], as.character)

factor_vars <- lapply(colData(RGset_controls), class) == "factor"
colData(RGset_controls)[, factor_vars] <- lapply(colData(RGset_controls)[, factor_vars], as.character)

# Juntando os dados dos arrays
rgset_all <- combineArrays(RGset_controls,RGset_HRC,outType = ("IlluminaHumanMethylation450k"))

rgset_all
# class: RGChannelSetExtended 
# dim: 575385 188 
# rownames(575385): 10600322 10600328 ... 74810485 74810492
# colnames(188): 1 2 ... I236_W1 I378_W1
# colData names(60): Methylation.Chip.ID CAG_ID ... AGE ArrayTypes
# Annotation
# array: IlluminaHumanMethylation450k
# annotation: ilmn12.hg19

############################################################
#################### QUALITY CONTROL  ######################
############################################################
# calculate the detection p-values
detP <- detectionP(rgset_all)
head(detP[,1:5])

# Importando informacoes fenotipicas participantes HRC + controles
target <- read.csv2(file = "~/Documents/ILS_INPD_48/Comparacao_Controles_CAG/targets.csv", header = T)
str(target)
target$sampleID <- as.character(target$sampleID)
a <- target$sampleID
b <- colnames(rgset_all)
identical(a,b) # [1] TRUE
gender <- as.numeric(target$GENDER)

# Examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(target$Group)], las=2,
        cex.names=0.8,ylab="Mean detection p-values")
abline(h=0.01,col="red")
legend("topleft", legend=levels(factor(target$Group)), fill=pal,
       bg="white")
barplot(colMeans(detP), col=pal[factor(target$Group)], las=2,
        cex.names=0.8, ylim = c(0,0.002), ylab="Mean detection p-values")
legend("topleft", legend=levels(factor(target$Group)), fill=pal,
       bg="white")

qcReport(rgset_all, sampNames=target$sampleID, sampGroups=target$Group,
         pdf="qcReport.pdf")

# remove poor quality samples
keep <- colMeans(detP) < 0.05
rgset_all <- rgset_all[,keep]
rgset_all
# class: RGChannelSetExtended 
# dim: 575385 188 

############################################################
#################### NORMALISATION  ########################
############################################################
# Normalization by preprocessFunnorm (faz tbm a NOOB background correction)
rgset_norm <- preprocessFunnorm(rgset_all, nPCs=2, sex = gender, bgCorr = TRUE,
                                dyeCorr = TRUE, keepCN = TRUE, ratioConvert = TRUE, 
                                verbose = TRUE)

# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(rgset_all, sampGroups=target$Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(target$Group)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(rgset_norm), sampGroups=target$Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(target$Group)),
       text.col=brewer.pal(8,"Dark2"))

############################################################
################### DATA EXPLORATION #######################
############################################################
# MDS plots to look at largest sources of variation
pal <- brewer.pal(8,"Dark2")
plotMDS(getM(rgset_norm), top=1000, gene.selection="common",
        col=pal[factor(target$Group)])
legend("top", legend=levels(factor(target$Group)), text.col=pal,
       bg="white", cex=0.7)
plotMDS(getM(rgset_norm), top=1000, gene.selection="common",
        col=pal[factor(target$GENDER)])
legend("top", legend=levels(factor(target$GENDER)), text.col=pal,
       bg="white", cex=0.7)

# Examine higher dimensions to look at other sources of variation
plotMDS(getM(rgset_norm), top=1000, gene.selection="common",
        col=pal[factor(target$Group)], dim=c(1,3))
legend("top", legend=levels(factor(target$Group)), text.col=pal,
       bg="white", cex=0.7)
plotMDS(getM(rgset_norm), top=1000, gene.selection="common",
        col=pal[factor(target$Group)], dim=c(2,3))
legend("top", legend=levels(factor(target$Group)), text.col=pal,
       bg="white", cex=0.7)
plotMDS(getM(rgset_norm), top=1000, gene.selection="common",
        col=pal[factor(target$Group)], dim=c(3,4))
legend("top", legend=levels(factor(target$Group)), text.col=pal,
       bg="white", cex=0.7)

############################################################
###################### FILTERING  ##########################
############################################################
# ensure probes are in the same order in the rgset_norm and detP objects
detP <- detP[match(featureNames(rgset_norm),rownames(detP)),]

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(rgset_norm)
table(keep)
# FALSE   TRUE 
# 8477  444355 

rgset_norm_Flt <- rgset_norm[keep,]
rgset_norm_Flt
# class: GenomicRatioSet
# dim: 444355 188 

# if your data includes males and females, remove probes on the sex chromosomes
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

keep <- !(featureNames(rgset_norm_Flt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])

table(keep)
# FALSE   TRUE 
# 9909    434446 

rgset_norm_Flt <- rgset_norm_Flt[keep,]
rgset_norm_Flt
# class: GenomicRatioSet 
# dim: 434446 188 

# remove probes with SNPs at CpG site
rgset_norm_Flt <- dropLociWithSnps(rgset_norm_Flt)
rgset_norm_Flt
# class: GenomicRatioSet 
# dim: 420347 188 

# exclude cross reactive probes
# OBSERVACAO: non-CpG- targeting probes (Probe ID prefix = “ch”) 
xReactiveProbes <- read.csv2(file=paste(dataDirectory,
                                        "Comparacao_Controles_CAG/48639-non-specific-probes-Illumina450k.csv",
                                        sep="/"))
keep <- !(featureNames(rgset_norm_Flt) %in% xReactiveProbes$TargetID)
table(keep)
# FALSE   TRUE 
# 24,495 395,852

rgset_norm_Flt <- rgset_norm_Flt[keep, ]
rgset_norm_Flt
# class: GenomicRatioSet 
# dim: 395852 188

# calculate M-values for statistical analysis
mvalues <- getM(rgset_norm_Flt)

# remove probes with SNPs at CpG site according to DMRcate package
mvalues <- rmSNPandCH(mvalues, dist=2, mafcut=0.1, rmXY = TRUE, rmcrosshyb = TRUE)
dim(mvalues) # [1] 389,662    188

# Excluir amostras 36, 55 (outlier nos graficos) 94
# 36: report = female; no grafico se agrupa como male
mvalues <- mvalues[,-c(36,55,94)]
dim(mvalues) # [1] 389,662    185

target <- target[-c(36,55,94), ]
dim(target) # [1] 185  10

a <- target$sampleID
b <- colnames(mvalues)
identical(a,b) # [1] TRUE

############################################################
################ EWAS USING AGE as IV ######################
############################################################
# IV: independent variable

# Selecionando so os controles e excluindo amostras do HRC 
mvalues_age <- mvalues[,1:137]
target_age <- target[1:137, ]

a <- target_age$sampleID
b <- colnames(mvalues_age)
identical(a,b) # [1] TRUE

# Variaveis no modelo
age <- target_age$AGE

# Modelo estatistico
design <- model.matrix(~age)
head(design)

myannotation <- cpg.annotate(datatype = c("array"), mvalues_age, what=c("M"),arraytype=c("450K"), 
                             analysis.type="differential", design = design, fdr = 0.05, coef="age")
# Your contrast returned 646 individually significant probes

cpgs <- data.frame(ID = myannotation$ID, weights = abs(myannotation$stat),CHR = as.character(myannotation$CHR),
                   pos = myannotation$pos,betafc = myannotation$betafc,indfdr = myannotation$indfdr,
                   is.sig = myannotation$is.sig)

# Anotacoes genomicas do 450K
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- as.data.frame(ann450k)
dim(ann450k) # [1] 485,512     33
dim(cpgs) # [1] 389,662      7

cpgs$ID <- as.character(cpgs$ID)
ann450k$ID <- rownames(ann450k)

cpgs_anotado <- left_join(cpgs,ann450k, by = "ID")
dim(cpgs_anotado) # [1] 389,662     40

# Salvando cpgs_anotado:
write.csv2(cpgs_anotado, file = "Associacao_Idade_TODAS_Informacoes.csv", 
           quote = F, row.names = F)
write.table(cpgs_anotado,file="Associacao_Idade_TODAS_Informacoes.txt", sep="\t", 
            row.names=F,quote=FALSE)

############################################################
######## Excluding Probes Age and Puberty Associated  ######
############################################################
# Probes to remove: Todas as probes associadas com idade
# na analise anterior utilizando os controles do Hakon +
# Probes encontradas pelo Almstrup et al., 2016 como

# Importando os IDs das probes que serao removidas
probes_to_remove <- read.csv2(file = "CpG_Ids_ParaRemover.csv",header = T)
probes_to_remove$ID <- as.character(probes_to_remove$ID)
dim(probes_to_remove) # [1] 757   1

# Selecionando so as informacoes do HRC
mvalues_b <- mvalues[,138:185]
dim(mvalues_b) # [1] 389,662    48

target_b <- target[138:185, ]
dim(target_b) # [1] 48  10

a <- target_b$sampleID
b <- colnames(mvalues_b)
identical(a,b) # [1] TRUE

# Removendo de Mvalues as probes associadas com idade e puberdade
remove <- probes_to_remove$ID

mvalues_b <- as.data.frame(mvalues_b)
mvalues_b$ID <- rownames(mvalues_b)

mvalues_b <- mvalues_b[!(mvalues_b$ID %in% remove), ]
dim(mvalues_b) # [1] 388,924     49
# 389,662 - 388,924 = 738 probes

mvalues_b <- mvalues_b[,1:48]
mvalues_b <- as.matrix(mvalues_b)

############################################################
################# Finding DMP and DMR ######################
############################################################
# Analise so com os dados do HRC 

# Como targets_HRC e target_b sao iguais, da para usar as informacoes fenotipicas das duas data frames
a <- as.character(targets_HRC$sampleID)
b <- target_b$sampleID
identical(a,b) # TRUE

# Coeficients
time <- factor(targets_HRC$TIME, levels=c("Baseline", "Followup"))
pair <- factor(targets_HRC$pair)
# Covariables
chip_number <- targets_HRC$CHIP_NUMBER
chip_position <- as.numeric(targets_HRC$CHIP_POSITION)

# Modelo
design <- model.matrix(~time+pair+chip_number+chip_position)
head(design)

######## DMP ######
myannotation <- cpg.annotate(datatype = c("array"), mvalues_b, what=c("M"),arraytype=c("450K"), analysis.type="differential", 
                             design = design, fdr = 0.05, coef="timeFollowup")
# Your contrast returned 663 individually significant probes.

######## DMR ######
DMRoutput <- dmrcate(myannotation, lambda=1000, C=2,p.adjust.method = "BH", min.cpgs = 3)
dim(DMRoutput$results) # [1] 90  6
DMR <- DMRoutput$results

write.csv2(DMR, file = "DMR_Longitudinal.csv", 
           quote = F, row.names = F)
write.table(DMR,file="DMR_Longitudinal.txt", sep="\t", 
            row.names=F,quote=FALSE)

############################################################
################ GENE EXPRESSION ANALYSIS  #################
############################################################
# Importando os dados de expressao genica
load("Dados_Expressao_eset_final.RData")
eset_final # 6782 features, 228 samples 

# selecionando os dados do grupo incidente
pdat <- pData(eset_final)
pdat <- pdat[order(pdat$Sample_Group), ]
pdat_inc <- pdat[58:109,]

# 20214,20228 e 21044 esta a mais
# tirar 20795 e 20815
pdat_inc <- pdat_inc[order(pdat_inc$Sample_ID), ]
pdat_met <- pdat_inc[-c(9:12,49,50,33,34), ]

# Selecionando os dados das probes
met <- as.character(pdat_met$sampleID)
eset_met <- eset_final[,met]
dim(eset_met)
# Features  Samples 
# 6782       44 

probes <- fData(eset_final)
genes_no_array <- probes[,c(1:3)]
write.csv2(genes_no_array, file = "genes_que_tem_no_GENEX.csv")

# CpGs
genes <- read.csv2(file = "Genes_DifMet_ComExpressao_GENEX.csv", header = T)
dim(genes) # [1] 122   3
length(unique(genes$UCSC_RefGene_Name)) # [1] 103

# Usando table para verificar se tem duplicacao:
n_occur <- data.frame(table(genes$UCSC_RefGene_Name)) # me da uma tabela com a lista de IDs e o numero de vezes que elas ocorrem
n_occur[n_occur$Freq > 1,] # mostra quais sao os IDs duplicados
genes[genes$UCSC_RefGene_Name %in% n_occur$Var1[n_occur$Freq > 1],] # mostra as linhas duplicadas da data frame

genes_row <- as.character(genes$Tem_Expressao_GENEX_Rowname) # 122 caracteres

# Dados brutos de expressao genica
d_matrix <- exprs(eset_met)
d_matrix <- as.data.frame(d_matrix)

# selecionando os dados dos genes diferencialmente metilados do artigo 1
d_matrix_met <- d_matrix[genes_row, ]
d_matrix_met <- as.matrix(d_matrix_met)
dim(d_matrix_met) # [1] 122  44

group <- factor(pdat_met$PHENOTYPE, levels=c("A", "D"))
levels(group)[levels(group)=="A"] <- "Baseline"
levels(group)[levels(group)=="D"] <- "Follow-up"

n <- 2
pair <- as.factor(rep(1:22,each=n))

# Modelo Estatistico
design <- model.matrix(~group+pair) 
head(design)

fit <- lmFit(d_matrix_met, design)
fit <- eBayes(fit)

results <- topTable(fit, coef="groupFollow-up", adjust.method="BH", number = nrow(d_matrix_met))

dim(results) # [1] 122   7
results$probe_row <- rownames(results)
genes$UCSC_RefGene_Name <- as.character(genes$UCSC_RefGene_Name)
genes$Tem_Expressao_GENEX_Rowname <- as.character(genes$Tem_Expressao_GENEX_Rowname)
names(genes)[2] <- "probe_row"

results <- full_join(results,genes,by="probe_row")
results <- results[,c(1:6,8,9,7)]

# Salvando cpg_genes:
write.csv2(results, file = "GENEX_longitudunal.csv", 
           quote = F, row.names = T)
write.table(results,file="GENEX_longitudunal.txt", sep="\t", 
            row.names=T,quote=FALSE)