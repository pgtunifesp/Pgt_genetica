########################################################################
# Title: Gene Expression, DNA Methylation and Environment
# Author: Leticia Spindola
# Date: 23/05/2018
# Purpose: Analise dos dados de metilacao, expressao e dados ambientais
#######################################################################

library(minfi)
require(dplyr)
library(ggplot2)
library(limma)

dataDirectory <- "~/Documents/ILS_INPD_48"
setwd("~/Documents/ILS_INPD_48")

############################################################
################### Dados de Expressao  ####################
############################################################
# Importando os dados de expressao genica
eset_met <- readRDS("eset_met.RDS")

# salvando os dados brutos de expressao em uma matrix
d_matrix <- exprs(eset_met)
pdata <- pData(eset_met)

# arrumar_IDs <- colnames(d_matrix)
# sampleID <- pdata$sampleID
# identical(arrumar_IDs,sampleID) # [1] TRUE

colnames(d_matrix) <- c("I211_W0","I211_W1","I175_W0","I175_W1","I070_W0","I070_W1","I051_W0",
                        "I051_W1","I117_W0","I117_W1","I144_W0","I144_W1","I069_W0","I069_W1",
                        "I197_W0","I197_W1","I209_W0","I209_W1","I227_W0","I227_W1","I224_W0",
                        "I224_W1","I242_W0","I242_W1","I310_W0","I310_W1","I236_W0","I236_W1",
                        "I378_W0","I378_W1","I272_W0","I272_W1","I281_W0","I281_W1","I338_W0",
                        "I338_W1","I339_W0","I339_W1","I309_W0","I309_W1","I330_W0","I330_W1",
                        "I331_W0","I331_W1")

probes <- fData(eset_met)

# Importando a lista para o R
DEG <- read.csv2(file = "Genes_diff_expressos_Incidentes2.csv",header = F)
# View(DEG)
# class(DEG) # [1] "data.frame"

# Transformando em um vetor
DEG <- as.character(DEG$V1)
DEG

# Tirando genes duplicados
DEG[duplicated(DEG)] # character(0)
length(DEG) # [1] 66

# Transformando as informacoes dos genes em uma string so
g <- paste(DEG,collapse = "|") 

# Criando uma coluna com informacao de row.names em probes 
probes$row_names <- rownames(probes)

# Selecionando so as probes dos genes diferencialmente expressos nos incidentes
genes_exp_incidentes <- filter(probes,grepl(g, ProbeID))
dim(genes_exp_incidentes) # [1] 72 30
x <- genes_exp_incidentes$SYMBOL

x[duplicated(x)] 
# [1] "PEX16"

keep <- genes_exp_incidentes$row_names

d_matrix_select <- d_matrix[keep,]
dim(d_matrix_select) # [1] 72 44

### Vector necessario para as analises: identificando quem eh quem (pair)
n <- 2
pair<- as.factor(rep(1:22,each=n))

############################################################
################### Dados de Metilacao  ####################
############################################################
# Importando os dados do EPIC 
rgset_HRC_norm_Flt <- readRDS("rgset_HRC_norm_Flt.RDS")

# Pegando as informacoes de annotacao de probe do EPIC
annotation <- rgset_HRC_norm_Flt %>% getAnnotation %>% as.data.frame

# Importando as informacoes de phenotype
targets_HRC <- read.csv2(file = "/Users/leticiaspindola/Documents/ILS_INPD_48/pdata_INPDmeth_final_25_07_2018.csv", 
                         header = T)
targets_HRC <- targets_HRC[order(targets_HRC$sampleID),] # # ordenando as linhas do genetic ID
targets_HRC_44_total <- targets_HRC[-c(27,28,33,34),] # remover as linhas 27 e 28 (I250) e 33 e 34 (I295)

rownames(targets_HRC_44_total) <- targets_HRC_44_total$sampleID
targets_HRC_44 <-targets_HRC_44_total [,c("Un_ev_2_NEW","Inter_2_NEW","Cont_chan_2_NEW",
                                          "School_2_NEW","Health_2_NEW","Env_stress_2_NEW","TIME")]
targets_HRC_44$GeneticID <- rownames(targets_HRC_44)
targets_HRC_44 <- filter(targets_HRC_44,TIME=="Followup")
targets_HRC_44 <- select(targets_HRC_44,-TIME)
rownames(targets_HRC_44) <- targets_HRC_44$GeneticID
targets_HRC_44 <- select(targets_HRC_44,-GeneticID)

# calculate M-values for statistical analysis
mvalues <- getM(rgset_HRC_norm_Flt)
mvalues <- mvalues[,order(colnames(mvalues))] # ordenando as colunas dos dados de metilacao
mvalues <- mvalues[,-c(27,28,33,34)] # remover as linhas 27 e 28 (I250) e 33 e 34 (I295)

# calculate Beta-values for plot and interpretation of statistical analysis
betas <- getBeta(rgset_HRC_norm_Flt)
betas <- betas[,order(colnames(betas))] # ordenando as colunas dos dados de metilacao
betas <- betas[,-c(27,28,33,34)] # remover as linhas 27 e 28 (I250) e 33 e 34 (I295)

# Ordenado as colunas dos dados de expressao
d_matrix_select <- d_matrix_select[,order(colnames(d_matrix_select))]

a <- colnames(mvalues)
b <- colnames(d_matrix_select)
c <- colnames(betas)
d <- as.character(targets_HRC_44_total$sampleID)
identical(a,b) # [1] TRUE
identical(a,c) # [1] TRUE
identical(b,c) # [1] TRUE
identical(a,d) # [1] TRUE
identical(b,d) # [1] TRUE
identical(c,d) # [1] TRUE

# Para as analises abaixo, preciso dos objetos:
# annotation, mvalues, betas, d_matrix_select, targets_HRC_44 e targets_HRC_44_total
rm(d_matrix,eset_met,pdata,probes,rgset_HRC_norm_Flt,targets_HRC)
gc()

############ Selecionando CpGs do gene SERTAD2 #############
# Verificando quantas CpGs tem no gene analisado
SERTAD2_cpgs <- filter(annotation,grepl("NM_014755", UCSC_RefGene_Accession))
dim(SERTAD2_cpgs) # [1] 28 46

keep <- SERTAD2_cpgs$Name

SERTAD2_mvalues <- mvalues[keep,]
dim(SERTAD2_mvalues) # [1] 28 44

# Pegando os dados de expressao genica
SERTAD2 <- as.numeric(d_matrix_select["xorMF.VUSF4ie5f9zk",])

# y = b + ax
# y = dados de expressao
# x = dados de metilacao

# Criando FOR
# Criar uma data frame com 6 colunas e nrow(SERTAD2_mvalues) linhas:
n <- nrow(SERTAD2_mvalues)
SERTAD2_exp_met <- data.frame(id=rep(0, n), estimate=rep(0, n), SE=rep(0, n),
                              t=rep(0, n), p=rep(0, n), r.sq=rep(0, n))

for (i in 1:nrow(SERTAD2_mvalues)) { # mudar aqui
  row <- SERTAD2_mvalues[i,] # mudar aqui
  id <- rownames(SERTAD2_mvalues)[i] # mudar aqui
  x <- lm(SERTAD2~row+pair)  # Mudar aqui
  estimate <- summary(x)$coefficients[2,1]
  SE <- summary(x)$coefficients[2,2]
  t <- summary(x)$coefficients[2,3]
  p <- summary(x)$coefficients[2,4]
  r.sq <- summary(x)$r.squared
  SERTAD2_exp_met[i,] <- c(id,estimate,SE,t,p,r.sq) # mudar aqui
}

# Trasformar as colunas estimate, SE, t, p e r.sq de character para numeric
cols <- c(2, 3, 4, 5, 6)    
SERTAD2_exp_met[,cols] <- apply(SERTAD2_exp_met[,cols], 2, function(x) as.numeric(as.character(x)))

# Correcao de BONFERRONI
SERTAD2_exp_met <- SERTAD2_exp_met[order(SERTAD2_exp_met$p),]
SERTAD2_exp_met$adj.p <- SERTAD2_exp_met$p * nrow(SERTAD2_exp_met)

# Salvando cpg_genes:
write.csv2(SERTAD2_exp_met, file = "Paper2_SERTAD2_y_EXP_x_METH_mvalues.csv", 
           quote = F, row.names = T)
write.table(SERTAD2_exp_met,file="Paper2_SERTAD2_y_EXP_x_METH_mvalues.txt", sep="\t", 
            row.names=T,quote=FALSE)

############ Metilacao e Trauma ####
# Verificando se metilacao de cg10184289 foi influenciada por trauma
# y = b + ax
# y = dados de metilacao (delta)
# x = dados de trauma
SERTAD2_cg <- as.data.frame(SERTAD2_mvalues["cg10184289",])
colnames(SERTAD2_cg) <- c("SERTAD2_cg")
SERTAD2_cg$pair <- pair 
SERTAD2_cg$ID <- rownames(SERTAD2_cg)

SERTAD2_cg <- SERTAD2_cg %>% 
  group_by(pair) %>% 
  mutate(delta_cg = SERTAD2_cg - lag(SERTAD2_cg))

SERTAD2_cg <- SERTAD2_cg[,c(3,4)]
SERTAD2_cg <- na.omit(SERTAD2_cg)
rownames(SERTAD2_cg) <- SERTAD2_cg$ID

a <- rownames(SERTAD2_cg)
b <- rownames(targets_HRC_44)
identical(a,b) # [1] TRUE

SERTAD2_cg <- as.character(SERTAD2_cg$delta_cg)

n <- ncol(targets_HRC_44)
SERTAD2_meth_trauma <- data.frame(id=rep(0, n), 
                                  estimate=rep(0, n),
                                  SE=rep(0, n),
                                  t=rep(0, n),
                                  p=rep(0, n),
                                  r.sq=rep(0, n))

for (i in 1:ncol(targets_HRC_44)) {
  col <- targets_HRC_44[,i]
  id <- colnames(targets_HRC_44)[i]
  x <- lm(SERTAD2_cg~col) # mudar aqui
  estimate <- summary(x)$coefficients[2,1]
  SE <- summary(x)$coefficients[2,2]
  t <- summary(x)$coefficients[2,3]
  p <- summary(x)$coefficients[2,4]
  r.sq <- summary(x)$r.squared
  SERTAD2_meth_trauma[i,] <- c(id,estimate,SE,t,p,r.sq)
}

# Trasformar as colunas estimate, SE, t, p e r.sq de character para numeric
SERTAD2_meth_trauma[,cols] = apply(SERTAD2_meth_trauma[,cols], 2, function(x) as.numeric(as.character(x)))

# Correcao de BONFERRONI
SERTAD2_meth_trauma <- SERTAD2_meth_trauma[order(SERTAD2_meth_trauma$p),]
SERTAD2_meth_trauma$adj.p <- SERTAD2_meth_trauma$p * nrow(SERTAD2_meth_trauma)
SERTAD2_meth_trauma

# Salvando:
write.csv2(SERTAD2_meth_trauma, file = "Paper2_SERTAD2_y_METH_x_TRAUMA_mvalues.csv", 
           quote = F, row.names = T)
write.table(SERTAD2_meth_trauma,file="Paper2_SERTAD2_y_METH_x_TRAUMA_mvalues.txt", sep="\t", 
            row.names=T,quote=FALSE)

rm(SERTAD2,SERTAD2_exp_met,SERTAD2_cpgs,
   SERTAD2_meth_trauma,SERTAD2_mvalues,SERTAD2_cg10184289,SERTAD2_cg)
gc()
