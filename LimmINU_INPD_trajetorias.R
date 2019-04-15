#Bibliotecas necessarias
library(oligo)
library(limma)
library(bioDist)
library(annotate)
library(hugene10sttranscriptcluster.db)
library(lumi)
library(xtable)
library(plyr)
library(CellMix)

###########
# Análise para identificar probes diferencialmente expressas entre Caso x Controle (CACO) 
# e entre Follow-up x Casos (FOCA)
###########

setwd("~/Documents/GENEX_INPD/INPD_1st_lumi_processing_03_02_2017_1486127173/eset_final/")

#Chamar os resultado do pipleine "foca2_corrigido.R" (adaptado de S.NEWHOUSE no github)
load("INPD_1st.eset_final.RData")
eset_final

# Banco atualizado com as caractrísticas da amostra 
pDat <- eset_final@phenoData@data

# dados clínicos
INPDw2 <- read.csv("~/Documents/GENEX_INPD/HRC_DB_reportnumber_3_080517forMarcos.csv",
header=TRUE, stringsAsFactors=FALSE, fileEncoding="latin1")
names(INPDw2)


INPDw2$Sample_ID <- INPDw2$ï..subjectid

#dados ged sangue

teste3 <- gedBlood(eset_final)

dados<- .coef(teste3)

blood <- as.data.frame(t(dados))
blood$Sample_ID <- row.names(blood)

# juntar todos
juntado <- join_all(list(pDat,INPDw2), by= "Sample_ID", match = "all")

eset_final@phenoData@data <- juntado

# Grupos de comparacao. 
# Nao longitudinal > Co = Controls; Ca = Cases, Incidentes = In, Remitentes = Re
# FIXED = ser aquele fenótipo nos wave 0 e wave1
# state = ser aquele fenótipo naquele momento, ex.: quem está no grupo incidente no wave 0 era "controle" e no wave 1 caso

###FIXED CoCa 0w: analise caso x controle no wave 0
Coca_0w <- pDat[pDat$Sample_Group=="A"|pDat$Sample_Group=="D",]
caquito_0w <- as.character(Coca_0w$sampleID)
eset_caco_0w <- eset_final[,caquito_0w]
dim(eset_caco_0w)

###FIXED CoCa 1w: analise caso x controle no wave 1
Coca_1w <- pDat[pDat$Sample_Group=="Ap"|pDat$Sample_Group=="Dp",]
caquito_1w <- as.character(Coca_1w$sampleID)
eset_caco_1w <- eset_final[,caquito_1w]
dim(eset_caco_1w)

###Rein 0w: analise remitente x incidente no wave 0
Coca_0w <- pDat[pDat$Sample_Group=="B"|pDat$Sample_Group=="C",]
caquito_0w <- as.character(Coca_0w$sampleID)
eset_rein_0w <- eset_final[,caquito_0w]
dim(eset_rein_0w)

###Rein 1w: analise remitente x incidente no wave 1
Coca_1w <- pDat[pDat$Sample_Group=="Bp"|pDat$Sample_Group=="Cp",]
caquito_1w <- as.character(Coca_1w$sampleID)
eset_Rein_1w <- eset_final[,caquito_1w]
dim(eset_Rein_1w)

#####CoCa state 0W
Coca_0w <- pDat[pDat$Sample_Group=="A"|pDat$Sample_Group=="C"|pDat$Sample_Group=="B"|pDat$Sample_Group=="D",]
caquito_0w <- as.character(Coca_0w$sampleID)
eset_CocaST_0w <- eset_final[,caquito_0w]
dim(eset_CocaST_0w)

###CoCa state 1w
Coca_1w <- pDat[pDat$Sample_Group=="Ap"|pDat$Sample_Group=="Cp"|pDat$Sample_Group=="Bp"|pDat$Sample_Group=="Dp",]
caquito_1w <- as.character(Coca_1w$sampleID)
eset_CocaST_1w <- eset_final[,caquito_1w]
dim(eset_CocaST_1w)

###Any Disorder = DAWBA desconsiderando a nossa selecao pelo CBCL
AnyEControls <- juntado[(juntado$Sample_Group=="Ap"&juntado$dcFUPany==0&juntado$GROUPS==2)| (juntado$dcFUPany==2&juntado$GROUPS==2),]
AnyEControls_1w <- as.character(AnyEControls$sampleID)
eset_ANYDIS_1w <- eset_final[,AnyEControls_1w]
dim(eset_ANYDIS_1w)

### Drugs: dorgas... todos nós utilizamos! Ver com a Leticia que raio de eh o "RiskFUPASUseever", acho que eh ever used para qualquer droga
juntado2 <- juntado[complete.cases(juntado$RiskFUPASUseever),]
Drug <- juntado2[(juntado2$RiskFUPASUseever==1|juntado2$RiskFUPASUseever==0)&juntado2$GROUPS==2,]
Drug_1w <- as.character(Drug$sampleID)
eset_drug_1w <- eset_final[,Drug_1w]
dim(eset_drug_1w)


########################
# PAREADO
########################
# NOVOS comparando trajetória controle e incidente
foca <- pDat[pDat$Sample_Group=="A"|pDat$Sample_Group=="Ap"|pDat$Sample_Group=="B"|pDat$Sample_Group=="Bp",]
eset_foca <- eset_final[,as.character(foca$sampleID)]
dim(eset_foca)
foco <- foca[duplicated(foca$Sample_ID),]
ficam <- as.character(foco$Sample_ID)
eset_PAREconinci <- eset_foca[,pData(eset_foca)$Sample_ID %in% ficam]

# NOVOS comparando trajetória CASO-CASO e remitente
# TODOS
foca <- pDat[pDat$Sample_Group=="D"|pDat$Sample_Group=="Dp"|pDat$Sample_Group=="C"|pDat$Sample_Group=="Cp"|pDat$Sample_Group=="A"|pDat$Sample_Group=="Ap"|pDat$Sample_Group=="B"|pDat$Sample_Group=="Bp",]
eset_foca <- eset_final[,as.character(foca$sampleID)]
dim(eset_foca)
foco <- foca[duplicated(foca$Sample_ID),]
ficam <- as.character(foco$Sample_ID)
eset_PAREcasremin <- eset_foca[,pData(eset_foca)$Sample_ID %in% ficam]


#curso normal
foca <- pDat[pDat$Sample_Group=="A"|pDat$Sample_Group=="Ap",]
eset_foca <- eset_final[,as.character(foca$sampleID)]
dim(eset_foca)
foco <- foca[duplicated(foca$Sample_ID),]
ficam <- as.character(foco$Sample_ID)
eset_PAREcontroles <- eset_foca[,pData(eset_foca)$Sample_ID %in% ficam]

#Incidentes
foca <- pDat[pDat$Sample_Group=="B"|pDat$Sample_Group=="Bp",]
eset_foca <- eset_final[,as.character(foca$sampleID)]
dim(eset_foca)
foco <- foca[duplicated(foca$Sample_ID),]
ficam <- as.character(foco$Sample_ID)
eset_PAREIncidentes <- eset_foca[,pData(eset_foca)$Sample_ID %in% ficam]

#Remientes
foca <- pDat[pDat$Sample_Group=="C"|pDat$Sample_Group=="Cp",]
eset_foca <- eset_final[,as.character(foca$sampleID)]
dim(eset_foca)
foco <- foca[duplicated(foca$Sample_ID),]
ficam <- as.character(foco$Sample_ID)
eset_PARERemitentes <- eset_foca[,pData(eset_foca)$Sample_ID %in% ficam]

#Casos ...
foca <- pDat[pDat$Sample_Group=="D"|pDat$Sample_Group=="Dp",]
eset_foca <- eset_final[,as.character(foca$sampleID)]
dim(eset_foca)
foco <- foca[duplicated(foca$Sample_ID),]
ficam <- as.character(foco$Sample_ID)
eset_PAREcasos <- eset_foca[,pData(eset_foca)$Sample_ID %in% ficam]

# baixo desvio
foca <- pDat[pDat$Sample_Group=="D"|pDat$Sample_Group=="Dp"|pDat$Sample_Group=="A"|pDat$Sample_Group=="Ap",]
eset_foca <- eset_final[,as.character(foca$sampleID)]
dim(eset_foca)
foco <- foca[duplicated(foca$Sample_ID),]
ficam <- as.character(foco$Sample_ID)
eset_PAREbaixodesvio <- eset_foca[,pData(eset_foca)$Sample_ID %in% ficam]

#Alto desvio
foca <- pDat[pDat$Sample_Group=="E"|pDat$Sample_Group=="Ep"|pDat$Sample_Group=="B"|pDat$Sample_Group=="Bc"|pDat$Sample_Group=="C"|pDat$Sample_Group=="Cp",]
eset_foca <- eset_final[,as.character(foca$sampleID)]
dim(eset_foca)
foco <- foca[duplicated(foca$Sample_ID),]
ficam <- as.character(foco$Sample_ID)
eset_PARE_altodesvio <- eset_foca[,pData(eset_foca)$Sample_ID %in% ficam]

#################################################################################
#################################################################################
#################################################################################
#################################################################################
# Variáveis que vou utilizar e que estão no meu pDat
grupos <- factor(phenoData(eset_caco_0w)$Sample_Group, levels=c("A", "D"))
sex <- factor(phenoData(eset_caco_0w)$SEX, levels=c("MALE", "FEMALE"))
# PC1 <- phenoData(eset_caco_0w)$PC1_all_pre_ComBat
B <- pDat$Ba # porcentagem dos tipos celulares para cada indivíduo (CellMix.R)
eri <- pDat$eri
mono <- pDat$monoa
#megak <- pDat$megak
NK2 <- pDat$NKa
gran <- pDat$grana
cd4 <- pDat$cd4a
cd8 <- pDat$cd8a
neutro <- pDat$NEUTRO
NK <- pDat$NK
monogb <- pDat$MONO
#cor <- pDat$COR
tc <- pDat$TC
tc_act <- pDat$TC_ACT
igg <- pDat$IGG
igm <- pDat$IGM
#dc <- pDat$DC_ACT
####################################
# quais serão as minhas covariaveis? 
####################################
#GedBlood
# Este foi utilizado o artigo 4
design <- model.matrix(~grupos+sex) 

# design <- model.matrix(~0+grupos)

#Fit
fit <- lmFit(eset_caco_0w, design)
fit2 <- contrasts.fit(fit, coef="gruposD")
fit3 <- eBayes(fit2)
topTable(fit3, coef="gruposD", adjust="bonferroni")
results <- decideTests(fit2)
results_bonferroni <- decideTests(fit2, adjust="bonferroni")

# Visualizar resultados com FDR e com bonferroni
#FDR
vennDiagram(results)
vennDiagram(results, include="up")
vennDiagram(results, include="down")
#Bonferroni
vennDiagram(results_bonferroni)
vennDiagram(results_bonferroni, include="up")
vennDiagram(results_bonferroni, include="down")

#hora de organizar os meus resultados
resultados_CACO <- topTable(fit3, coef="gruposD", adjust.method="bonferroni", 
                            number=Inf)

#volcano plot com o resultado de eBayes
volcanoplot(fit3, coef="gruposD", highlight=20, names=fit3$genes$SYMBOL)

#Só os top 100
#top1000_CACO <- topTable(fit3, coef="gruposCASE", adjust.method="bonferroni", number=1000)
top100_CACO <- topTable(fit3, coef="gruposD", adjust.method="bonferroni", number=100)

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:
write.table(resultados_CACO, file="Resultados_COCA0W.txt",sep="\t", 
            row.names=F,quote=FALSE)
write.table(top100_CACO,file="top100_COCA0W.txt", sep="\t", row.names=F,
            quote=FALSE)

########################
# PREPARAR PARA O WGCNA OU LEARNING MACHINE
########################
design <- model.matrix(~sex)
Fit1 <- lmFit(eset_caco, design)
correctedValues <- residuals.MArrayLM(Fit1, eset_caco)
write.table(correctedValues, file="Corrigido_MDS_cell_020216_GB.txt",sep=" ", row.names=F,
            quote=FALSE)

################################################################
#FIM DO COCA 0w
################################################################

#################################################################################
#################################################################################
#################################################################################
#################################################################################

################################################################
#INICIO DO COCA 1w
###############################################################
# Variáveis que vou utilizar e que estão no meu pDat
grupos <- factor(phenoData(eset_caco_1w)$Sample_Group, levels=c("Ap", "Dp")) #mudar
sex <- factor(phenoData(eset_caco_1w)$SEX, levels=c("MALE", "FEMALE")) #mudar
# PC1 <- phenoData(eset_caco_0w)$PC1_all_pre_ComBat

####################################
# quais serão as minhas covariaveis? 
####################################
#GedBlood
# Este foi utilizado o artigo 4
design <- model.matrix(~grupos+sex) 

# design <- model.matrix(~0+grupos)

#Fit
fit <- lmFit(eset_caco_1w, design) #mudar
fit2 <- contrasts.fit(fit, coef="gruposDp") #mudar
fit3 <- eBayes(fit2)
topTable(fit3, coef="gruposDp", adjust="bonferroni") #mudar
results <- decideTests(fit2)
results_bonferroni <- decideTests(fit2, adjust="bonferroni")

# Visualizar resultados com FDR e com bonferroni
#FDR
vennDiagram(results)
vennDiagram(results, include="up")
vennDiagram(results, include="down")
#Bonferroni
vennDiagram(results_bonferroni)
vennDiagram(results_bonferroni, include="up")
vennDiagram(results_bonferroni, include="down")

#hora de organizar os meus resultados
resultados_CACO <- topTable(fit3, coef="gruposDp", adjust.method="bonferroni", 
                            number=Inf) #mudar

#volcano plot com o resultado de eBayes
volcanoplot(fit3, coef="gruposDp", highlight=20, names=fit3$genes$SYMBOL) #mudar

#Só os top 100
#top1000_CACO <- topTable(fit3, coef="gruposCASE", adjust.method="bonferroni", number=1000)
top100_CACO <- topTable(fit3, coef="gruposDp", adjust.method="bonferroni", number=100) #mudar

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:
write.table(resultados_CACO, file="Resultados_COCA1W.txt",sep="\t", 
            row.names=F,quote=FALSE)
write.table(top100_CACO,file="top100_COCA1W.txt", sep="\t", row.names=F,
            quote=FALSE)

################################################################
#FIM DO COCA 1w
###############################################################

################################################################
#INICIO DO COCA ST 0W
###############################################################
# Variáveis que vou utilizar e que estão no meu pDat
grupos <- factor(phenoData(eset_CocaST_0w)$PHENOTYPE, levels=c("A", "D")) #mudar
sex <- factor(phenoData(eset_CocaST_0w)$SEX, levels=c("MALE", "FEMALE")) #mudar
# PC1 <- phenoData(eset_caco_0w)$PC1_all_pre_ComBat

####################################
# quais serão as minhas covariaveis? 
####################################
#GedBlood
# Este foi utilizado o artigo 4
design <- model.matrix(~grupos+sex) 

# design <- model.matrix(~0+grupos)

#Fit
fit <- lmFit(eset_CocaST_0w, design) #mudar
fit2 <- contrasts.fit(fit, coef="gruposD") #mudar
fit3 <- eBayes(fit2)
topTable(fit3, coef="gruposD", adjust="bonferroni") #mudar
results <- decideTests(fit2)
results_bonferroni <- decideTests(fit2, adjust="bonferroni")

# Visualizar resultados com FDR e com bonferroni
#FDR
vennDiagram(results)
vennDiagram(results, include="up")
vennDiagram(results, include="down")
#Bonferroni
vennDiagram(results_bonferroni)
vennDiagram(results_bonferroni, include="up")
vennDiagram(results_bonferroni, include="down")

#hora de organizar os meus resultados
resultados_CACO <- topTable(fit3, coef="gruposD", adjust.method="bonferroni", 
                            number=Inf) #mudar

#volcano plot com o resultado de eBayes
volcanoplot(fit3, coef="gruposD", highlight=20, names=fit3$genes$SYMBOL) #mudar

#Só os top 100
#top1000_CACO <- topTable(fit3, coef="gruposCASE", adjust.method="bonferroni", number=1000)
top100_CACO <- topTable(fit3, coef="gruposD", adjust.method="bonferroni", number=100) #mudar

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:
write.table(resultados_CACO, file="Resultados_COCAST0W.txt",sep="\t", 
            row.names=F,quote=FALSE)
write.table(top100_CACO,file="top100_COCAST0W.txt", sep="\t", row.names=F,
            quote=FALSE)

################################################################
#FIM DO COCA ST 0W
###############################################################

################################################################
#INICIO DO COCA ST 1W
###############################################################
# Variáveis que vou utilizar e que estão no meu pDat
grupos <- factor(phenoData(eset_CocaST_1w)$PHENOTYPE, levels=c("A", "D")) #mudar
sex <- factor(phenoData(eset_CocaST_1w)$SEX, levels=c("MALE", "FEMALE")) #mudar
# PC1 <- phenoData(eset_caco_0w)$PC1_all_pre_ComBat

####################################
# quais serão as minhas covariaveis? 
####################################
#GedBlood
# Este foi utilizado o artigo 4
design <- model.matrix(~grupos+sex) 

# design <- model.matrix(~0+grupos)

#Fit
fit <- lmFit(eset_CocaST_1w, design) #mudar
fit2 <- contrasts.fit(fit, coef="gruposD") #mudar
fit3 <- eBayes(fit2)
topTable(fit3, coef="gruposD", adjust="bonferroni") #mudar
results <- decideTests(fit2)
results_bonferroni <- decideTests(fit2, adjust="bonferroni")

# Visualizar resultados com FDR e com bonferroni
#FDR
vennDiagram(results)
vennDiagram(results, include="up")
vennDiagram(results, include="down")
#Bonferroni
vennDiagram(results_bonferroni)
vennDiagram(results_bonferroni, include="up")
vennDiagram(results_bonferroni, include="down")

#hora de organizar os meus resultados
resultados_CACO <- topTable(fit3, coef="gruposD", adjust.method="bonferroni", 
                            number=Inf) #mudar

#volcano plot com o resultado de eBayes
volcanoplot(fit3, coef="gruposD", highlight=20, names=fit3$genes$SYMBOL) #mudar

#Só os top 100
#top1000_CACO <- topTable(fit3, coef="gruposCASE", adjust.method="bonferroni", number=1000)
top100_CACO <- topTable(fit3, coef="gruposD", adjust.method="bonferroni", number=100) #mudar

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:
write.table(resultados_CACO, file="Resultados_COCAST1W.txt",sep="\t", 
            row.names=F,quote=FALSE)
write.table(top100_CACO,file="top100_COCAST1W.txt", sep="\t", row.names=F,
            quote=FALSE)

################################################################
#FIM DO COCA ST 1W
###############################################################

################################################################
#INICIO DO REIN 0W
###############################################################
# Variáveis que vou utilizar e que estão no meu pDat
grupos <- factor(phenoData(eset_rein_0w)$PHENOTYPE, levels=c("A", "D")) #mudar
sex <- factor(phenoData(eset_rein_0w)$SEX, levels=c("MALE", "FEMALE")) #mudar
# PC1 <- phenoData(eset_caco_0w)$PC1_all_pre_ComBat

####################################
# quais serão as minhas covariaveis? 
####################################
#GedBlood
# Este foi utilizado o artigo 4
design <- model.matrix(~grupos+sex) 

# design <- model.matrix(~0+grupos)

#Fit
fit <- lmFit(eset_rein_0w, design) #mudar
fit2 <- contrasts.fit(fit, coef="gruposD") #mudar
fit3 <- eBayes(fit2)
topTable(fit3, coef="gruposD", adjust="bonferroni") #mudar
results <- decideTests(fit2)
results_bonferroni <- decideTests(fit2, adjust="bonferroni")

# Visualizar resultados com FDR e com bonferroni
#FDR
vennDiagram(results)
vennDiagram(results, include="up")
vennDiagram(results, include="down")
#Bonferroni
vennDiagram(results_bonferroni)
vennDiagram(results_bonferroni, include="up")
vennDiagram(results_bonferroni, include="down")

#hora de organizar os meus resultados
resultados_CACO <- topTable(fit3, coef="gruposD", adjust.method="bonferroni", 
                            number=Inf) #mudar

#volcano plot com o resultado de eBayes
volcanoplot(fit3, coef="gruposD", highlight=20, names=fit3$genes$SYMBOL) #mudar

#Só os top 100
#top1000_CACO <- topTable(fit3, coef="gruposCASE", adjust.method="bonferroni", number=1000)
top100_CACO <- topTable(fit3, coef="gruposD", adjust.method="bonferroni", number=100) #mudar

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:
write.table(resultados_CACO, file="Resultados_rein0W.txt",sep="\t", 
            row.names=F,quote=FALSE)
write.table(top100_CACO,file="top100_rein0W.txt", sep="\t", row.names=F,
            quote=FALSE)

################################################################
#FIM DO rein0w
###############################################################
################################################################
#INICIO DO REIN 1W
###############################################################
# Variáveis que vou utilizar e que estão no meu pDat
grupos <- factor(phenoData(eset_Rein_1w)$PHENOTYPE, levels=c("A", "D")) #mudar
sex <- factor(phenoData(eset_Rein_1w)$SEX, levels=c("MALE", "FEMALE")) #mudar
# PC1 <- phenoData(eset_caco_0w)$PC1_all_pre_ComBat

####################################
# quais serão as minhas covariaveis? 
####################################
#GedBlood
# Este foi utilizado o artigo 4
design <- model.matrix(~grupos+sex) 

# design <- model.matrix(~0+grupos)

#Fit
fit <- lmFit(eset_Rein_1w, design) #mudar
fit2 <- contrasts.fit(fit, coef="gruposD") #mudar
fit3 <- eBayes(fit2)
topTable(fit3, coef="gruposD", adjust="bonferroni") #mudar
results <- decideTests(fit2)
results_bonferroni <- decideTests(fit2, adjust="bonferroni")

# Visualizar resultados com FDR e com bonferroni
#FDR
vennDiagram(results)
vennDiagram(results, include="up")
vennDiagram(results, include="down")
#Bonferroni
vennDiagram(results_bonferroni)
vennDiagram(results_bonferroni, include="up")
vennDiagram(results_bonferroni, include="down")

#hora de organizar os meus resultados
resultados_CACO <- topTable(fit3, coef="gruposD", adjust.method="bonferroni", 
                            number=Inf) #mudar

#volcano plot com o resultado de eBayes
volcanoplot(fit3, coef="gruposD", highlight=20, names=fit3$genes$SYMBOL) #mudar

#Só os top 100
#top1000_CACO <- topTable(fit3, coef="gruposCASE", adjust.method="bonferroni", number=1000)
top100_CACO <- topTable(fit3, coef="gruposD", adjust.method="bonferroni", number=100) #mudar

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:
write.table(resultados_CACO, file="Resultados_rein1W.txt",sep="\t", 
            row.names=F,quote=FALSE)
write.table(top100_CACO,file="top100_rein1W.txt", sep="\t", row.names=F,
            quote=FALSE)

################################################################
#FIM DO rein 1w
###############################################################

#################################################################################
#################################################################################
#################################################################################
#################################################################################

################################################################
#INICIO DOs Pareado Controles x incidentes
###############################################################
phenoData(eset_PAREconinci)$gruposs <- phenoData(eset_PAREconinci)$Sample_Group
phenoData(eset_PAREconinci)$gruposs <- gsub("p","",phenoData(eset_PAREconinci)$gruposs)

grupos <- factor(phenoData(eset_PAREconinci)$gruposs, levels=c("A", "B"))
PAREgrupos <- factor(phenoData(eset_PAREconinci)$GROUPS, levels=c("1", "2")) #mudar
ID_INPD <- factor(phenoData(eset_PAREconinci)$Sample_ID) #mudar

TS <- paste(grupos, PAREgrupos, sep=".")

TS <- factor(TS, levels=c("A.1","A.2","B.1","B.2"))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)

corfit <- duplicateCorrelation(eset_PAREconinci,design,block=phenoData(eset_PAREconinci)$Sample_ID)
corfit$consensus

fit <- lmFit(eset_PAREconinci,design,block=phenoData(eset_PAREconinci)$Sample_ID,correlation=corfit$consensus)

cm <- makeContrasts(
            W0xW1forA = A.2-A.1,
            W0xW1forB = B.2-B.1,
            W0xdoenca = B.1-A.1,
            W1xdoenca = B.2-A.2,
   levels=design)

fit2 <- contrasts.fit(fit,cm)
fit2 <- eBayes(fit2)

results <- decideTests(fit2)
vennDiagram(results)

teste <- as.data.frame(results@.Data)
table(teste)
teste$soma <- rowSums(teste)
teste2 <- teste[teste$soma==2|teste$soma==-2,]

alfa <- topTable(fit2, coef="W0xW1forB", number = Inf, adjust.method = "bonferroni")

alfalimpo <- alfa[!rownames(alfa)%in%rownames(teste2),]

write.table(alfa, file = "~/Documents/W0xW1forB_100718_AllGenes.txt",quote = F,row.names = F,col.names = T,sep = "\t")
write.table(alfalimpo, file = "~/Documents/W0xW1forB_100718_CleanGenes.txt",quote = F,row.names = F,col.names = T,sep = "\t")

################################################################
#INICIO DOs Pareado TODOS
################################################################
phenoData(eset_PAREcasremin)$gruposs <- phenoData(eset_PAREcasremin)$Sample_Group
phenoData(eset_PAREcasremin)$gruposs <- gsub("p","",phenoData(eset_PAREcasremin)$gruposs)

grupos <- factor(phenoData(eset_PAREcasremin)$gruposs, levels=c("A","B","C","D"))
PAREgrupos <- factor(phenoData(eset_PAREcasremin)$GROUPS, levels=c("1", "2")) #mudar
ID_INPD <- factor(phenoData(eset_PAREcasremin)$Sample_ID) #mudar

TS <- paste(grupos, PAREgrupos, sep=".")
TS

TS <- factor(TS, levels=c("A.1","A.2","B.1","B.2","C.1","C.2","D.1","D.2"))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)

corfit <- duplicateCorrelation(eset_PAREcasremin,design,block=phenoData(eset_PAREcasremin)$Sample_ID)
corfit$consensus

fit <- lmFit(eset_PAREcasremin,design,block=phenoData(eset_PAREcasremin)$Sample_ID,correlation=corfit$consensus)

cm <- makeContrasts(
  W0xW1forA = A.2-A.1,
  W0xW1forB = B.2-B.1,
  W0xW1forC = C.2-C.1,
  W0xW1forD = D.2-D.1,
  levels=design)

cm2 <makeContrasts(
  W0xCD = C.1-D.1,
  W1xCD = C.2-D.2,
  W0xCA = C.1-A.1,
  W1xCA = C.2-A.2,
  W0xDA = D.1-A.1,
  W1xDA = D.2-A.2,
  levels=design
)

fit2 <- contrasts.fit(fit,cm)
fit2 <- eBayes(fit2)

results <- decideTests(fit2)
vennDiagram(results)
# 
# teste <- as.data.frame(results@.Data)
# table(teste)
# teste$soma <- rowSums(teste)
# teste2 <- teste[teste$soma==2|teste$soma==-2,]

alfa <- topTable(fit2, coef="W0xW1forA", number = Inf)
beta <- topTable(fit2, coef="W0xW1forB", number = Inf)
gama <- topTable(fit2, coef="W0xW1forC", number = Inf)
delta <- topTable(fit2, coef="W0xW1forD", number = Inf)

# alfalimpo <- alfa[!rownames(alfa)%in%rownames(teste2),]

write.table(alfa, file = "~/Documents/W0xW1forA_100718_AllGenes.txt",quote = F,row.names = F,col.names = T,sep = "\t")
write.table(beta, file = "~/Documents/W0xW1forB_100718_AllGenes.txt",quote = F,row.names = F,col.names = T,sep = "\t")
write.table(gama, file = "~/Documents/W0xW1forC_100718_AllGenes.txt",quote = F,row.names = F,col.names = T,sep = "\t")
write.table(delta, file = "~/Documents/W0xW1forD_100718_AllGenes.txt",quote = F,row.names = F,col.names = T,sep = "\t")

# write.table(alfalimpo, file = "~/Documents/W0xW1forB_100718_CleanGenes.txt",quote = F,row.names = F,col.names = T,sep = "\t")



#GedBlood
# PAREdesign <- model.matrix(~ID_PEP+PAREgrupos+neutro+NK+monogb+igg+igm+tc+tc_act)

ID_INPD <- factor(phenoData(eset_PAREcontroles)$Sample_ID) #mudar
PAREgrupos <- factor(phenoData(eset_PAREcontroles)$GROUPS, levels=c("1", "2")) #mudar

PAREdesign <- model.matrix(~PAREgrupos+ID_INPD)


PAREfit <- lmFit(eset_PAREcontroles, PAREdesign) #mudar

PAREfit2 <- eBayes(PAREfit)
# PAREresults <- decideTests(PAREfit2)
# PAREresults2 <- decideTests(PAREfit2, adjust="bonferroni")

PAREresultados <- topTable(PAREfit2,coef="PAREgrupos2", adjust="bonferroni", 
                           number=Inf) #mudar
PAREresultados <- PAREresultados[,c("SYMBOL","TRANSCRIPT", "ACCESSION","REFSEQ_ID",
                                    "ENTREZ_GENE_ID", "PROBE_ID","CHROMOSOME",
                                    "PROBE_CHR_ORIENTATION",
                                    "logFC","P.Value","adj.P.Val","B")]

volcanoplot(PAREfit2, coef="PAREgruposAp", highlight=10, names=PAREfit2$genes$SYMBOL) #mudar

PAREtop100_FOCA <- topTable(PAREfit2, coef="PAREgruposAp", adjust="bonferroni", 
                            number=100) #mudar
PAREtop100_FOCA <- PAREtop100_FOCA[,c("SYMBOL","TRANSCRIPT", "ACCESSION","REFSEQ_ID",
                                      "ENTREZ_GENE_ID", "PROBE_ID","CHROMOSOME",
                                      "PROBE_CHR_ORIENTATION",
                                      "logFC","P.Value","adj.P.Val","B")]

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:

write.table(PAREtop100_FOCA,file="PAREtop100_Controles.txt", sep="\t", 
            row.names=F,quote=FALSE) #mudar
write.table(PAREresultados,file="PAREresultados_Controles.txt", sep="\t", 
            row.names=F,quote=FALSE) #mudar

################################################################
#INICIO DOs Pareado Controles ou curso normal
###############################################################
# identificar as duplas e o grupo de comparaç˜ãao

PAREgrupos <- factor(phenoData(eset_PAREcontroles)$Sample_Group, levels=c("A", "Ap")) #mudar
ID_INPD <- factor(phenoData(eset_PAREcontroles)$Sample_ID) #mudar

#GedBlood
# PAREdesign <- model.matrix(~ID_PEP+PAREgrupos+neutro+NK+monogb+igg+igm+tc+tc_act)
PAREdesign <- model.matrix(~ID_PEP+PAREgrupos)
PAREdesign <- model.matrix(~PAREgrupos+ID_INPD)

PAREfit <- lmFit(eset_PAREcontroles, PAREdesign) #mudar

PAREfit2 <- eBayes(PAREfit)
# PAREresults <- decideTests(PAREfit2)
# PAREresults2 <- decideTests(PAREfit2, adjust="bonferroni")

PAREresultados <- topTable(PAREfit2,coef="PAREgruposAp", adjust="bonferroni", 
                           number=Inf) #mudar
PAREresultados <- PAREresultados[,c("SYMBOL","TRANSCRIPT", "ACCESSION","REFSEQ_ID",
                                    "ENTREZ_GENE_ID", "PROBE_ID","CHROMOSOME",
                                    "PROBE_CHR_ORIENTATION",
                                    "logFC","P.Value","adj.P.Val","B")]

volcanoplot(PAREfit2, coef="PAREgruposAp", highlight=10, names=PAREfit2$genes$SYMBOL) #mudar

PAREtop100_FOCA <- topTable(PAREfit2, coef="PAREgruposAp", adjust="bonferroni", 
                            number=100) #mudar
PAREtop100_FOCA <- PAREtop100_FOCA[,c("SYMBOL","TRANSCRIPT", "ACCESSION","REFSEQ_ID",
                                      "ENTREZ_GENE_ID", "PROBE_ID","CHROMOSOME",
                                      "PROBE_CHR_ORIENTATION",
                                      "logFC","P.Value","adj.P.Val","B")]

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:

write.table(PAREtop100_FOCA,file="PAREtop100_Controles.txt", sep="\t", 
            row.names=F,quote=FALSE) #mudar
write.table(PAREresultados,file="PAREresultados_Controles.txt", sep="\t", 
            row.names=F,quote=FALSE) #mudar
################################################################
#FIM DOs Pareados Controles ou curso normal
###############################################################

################################################################
#INICIO DOs Pareado Incidentes
###############################################################
# identificar as duplas e o grupo de comparaç˜ãao

PAREgrupos <- factor(phenoData(eset_PAREIncidentes)$Sample_Group, levels=c("B", "Bp")) #mudar
ID_INPD <- factor(phenoData(eset_PAREIncidentes)$Sample_ID) #mudar

#GedBlood
# PAREdesign <- model.matrix(~ID_PEP+PAREgrupos+neutro+NK+monogb+igg+igm+tc+tc_act)
PAREdesign <- model.matrix(~ID_INPD+PAREgrupos)
PAREdesign <- model.matrix(~PAREgrupos+ID_INPD)

PAREfit <- lmFit(eset_PAREIncidentes, PAREdesign) #mudar

PAREfit2 <- eBayes(PAREfit)
# PAREresults <- decideTests(PAREfit2)
# PAREresults2 <- decideTests(PAREfit2, adjust="bonferroni")

PAREresultados <- topTable(PAREfit2,coef="PAREgruposBp", adjust="bonferroni", 
                           number=Inf) #mudar

topTable(PAREfit2,coef="PAREgruposBp", adjust="bonferroni", 
         number=Inf) #mudar

PAREresultados <- PAREresultados[,c("SYMBOL","TRANSCRIPT", "ACCESSION","REFSEQ_ID",
                                    "ENTREZ_GENE_ID", "PROBE_ID","CHROMOSOME",
                                    "PROBE_CHR_ORIENTATION",
                                    "logFC","P.Value","adj.P.Val","B")]

volcanoplot(PAREfit2, coef="PAREgruposBp", highlight=10, names=PAREfit2$genes$SYMBOL) #mudar

PAREtop100_FOCA <- topTable(PAREfit2, coef="PAREgruposBp", adjust="bonferroni", 
                            number=100) #mudar
PAREtop100_FOCA <- PAREtop100_FOCA[,c("SYMBOL","TRANSCRIPT", "ACCESSION","REFSEQ_ID",
                                      "ENTREZ_GENE_ID", "PROBE_ID","CHROMOSOME",
                                      "PROBE_CHR_ORIENTATION",
                                      "logFC","P.Value","adj.P.Val","B")]

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:

write.table(PAREtop100_FOCA,file="PAREtop100_Incidentes.txt", sep="\t", 
            row.names=F,quote=FALSE) #mudar
write.table(PAREresultados,file="PAREresultados_Incidentes.txt", sep="\t", 
            row.names=F,quote=FALSE) #mudar
################################################################
#FIM DOs Pareados Incidentes 
###############################################################

################################################################
#INICIO DOs Pareado Remitentes
###############################################################
# identificar as duplas e o grupo de comparaç˜ãao

PAREgrupos <- factor(phenoData(eset_PARERemitentes)$Sample_Group, levels=c("C", "Cp")) #mudar
ID_INPD <- factor(phenoData(eset_PARERemitentes)$Sample_ID) #mudar

#GedBlood
# PAREdesign <- model.matrix(~ID_PEP+PAREgrupos+neutro+NK+monogb+igg+igm+tc+tc_act)
PAREdesign <- model.matrix(~ID_INPD+PAREgrupos)
PAREdesign <- model.matrix(~PAREgrupos+ID_INPD)

PAREfit <- lmFit(eset_PARERemitentes, PAREdesign) #mudar

PAREfit2 <- eBayes(PAREfit)
# PAREresults <- decideTests(PAREfit2)
# PAREresults2 <- decideTests(PAREfit2, adjust="bonferroni")

PAREresultados <- topTable(PAREfit2,coef="PAREgruposCp", adjust="bonferroni", 
                           number=Inf) #mudar

topTable(PAREfit2,coef="PAREgruposCp", adjust="bonferroni", 
         number=Inf) #mudar

PAREresultados <- PAREresultados[,c("SYMBOL","TRANSCRIPT", "ACCESSION","REFSEQ_ID",
                                    "ENTREZ_GENE_ID", "PROBE_ID","CHROMOSOME",
                                    "PROBE_CHR_ORIENTATION",
                                    "logFC","P.Value","adj.P.Val","B")]

volcanoplot(PAREfit2, coef="PAREgruposCp", highlight=10, names=PAREfit2$genes$SYMBOL) #mudar

PAREtop100_FOCA <- topTable(PAREfit2, coef="PAREgruposCp", adjust="bonferroni", 
                            number=100) #mudar
PAREtop100_FOCA <- PAREtop100_FOCA[,c("SYMBOL","TRANSCRIPT", "ACCESSION","REFSEQ_ID",
                                      "ENTREZ_GENE_ID", "PROBE_ID","CHROMOSOME",
                                      "PROBE_CHR_ORIENTATION",
                                      "logFC","P.Value","adj.P.Val","B")]

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:

write.table(PAREtop100_FOCA,file="PAREtop100_Remitentes.txt", sep="\t", 
            row.names=F,quote=FALSE) #mudar
write.table(PAREresultados,file="PAREresultados_Remitentes.txt", sep="\t", 
            row.names=F,quote=FALSE) #mudar

################################################################
#FIM DOs Pareados Remitentes
###############################################################

################################################################
#INICIO DOs Pareado Casos
###############################################################
# identificar as duplas e o grupo de comparaç˜ãao

PAREgrupos <- factor(phenoData(eset_PAREcasos)$Sample_Group, levels=c("D", "Dp")) #mudar
ID_INPD <- factor(phenoData(eset_PAREcasos)$Sample_ID) #mudar

#GedBlood
# PAREdesign <- model.matrix(~ID_PEP+PAREgrupos+neutro+NK+monogb+igg+igm+tc+tc_act)
PAREdesign <- model.matrix(~ID_INPD+PAREgrupos)
PAREdesign <- model.matrix(~PAREgrupos+ID_INPD)

PAREfit <- lmFit(eset_PAREcasos, PAREdesign) #mudar

PAREfit2 <- eBayes(PAREfit)
# PAREresults <- decideTests(PAREfit2)
# PAREresults2 <- decideTests(PAREfit2, adjust="bonferroni")

PAREresultados <- topTable(PAREfit2,coef="PAREgruposDp", adjust="bonferroni", 
                           number=Inf) #mudar

topTable(PAREfit2,coef="PAREgruposDp", adjust="bonferroni", 
         number=Inf) #mudar

PAREresultados <- PAREresultados[,c("SYMBOL","TRANSCRIPT", "ACCESSION","REFSEQ_ID",
                                    "ENTREZ_GENE_ID", "PROBE_ID","CHROMOSOME",
                                    "PROBE_CHR_ORIENTATION",
                                    "logFC","P.Value","adj.P.Val","B")]

volcanoplot(PAREfit2, coef="PAREgruposDp", highlight=10, names=PAREfit2$genes$SYMBOL) #mudar

PAREtop100_FOCA <- topTable(PAREfit2, coef="PAREgruposDp", adjust="bonferroni", 
                            number=100) #mudar
PAREtop100_FOCA <- PAREtop100_FOCA[,c("SYMBOL","TRANSCRIPT", "ACCESSION","REFSEQ_ID",
                                      "ENTREZ_GENE_ID", "PROBE_ID","CHROMOSOME",
                                      "PROBE_CHR_ORIENTATION",
                                      "logFC","P.Value","adj.P.Val","B")]

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:

write.table(PAREtop100_FOCA,file="PAREtop100_Casos.txt", sep="\t", 
            row.names=F,quote=FALSE) #mudar
write.table(PAREresultados,file="PAREresultados_Casos.txt", sep="\t", 
            row.names=F,quote=FALSE) #mudar

################################################################
#FIM DOs Pareados Casos
###############################################################

################################################################
#INICIO DOs Pareado baixo desvio
###############################################################
# identificar as duplas e o grupo de comparaç˜ãao

PAREgrupos <- factor(phenoData(eset_PAREbaixodesvio)$GROUPS, levels=c("1", "2")) #mudar
ID_INPD <- factor(phenoData(eset_PAREbaixodesvio)$Sample_ID) #mudar

#GedBlood
# PAREdesign <- model.matrix(~ID_PEP+PAREgrupos+neutro+NK+monogb+igg+igm+tc+tc_act)
PAREdesign <- model.matrix(~ID_INPD+PAREgrupos)
PAREdesign <- model.matrix(~PAREgrupos+ID_INPD)

PAREfit <- lmFit(eset_PAREbaixodesvio, PAREdesign) #mudar

PAREfit2 <- eBayes(PAREfit)
# PAREresults <- decideTests(PAREfit2)
# PAREresults2 <- decideTests(PAREfit2, adjust="bonferroni")

PAREresultados <- topTable(PAREfit2,coef="PAREgrupos2", adjust="bonferroni", 
                           number=Inf) #mudar

topTable(PAREfit2,coef="PAREgrupos2", adjust="bonferroni", 
         number=Inf) #mudar

PAREresultados <- PAREresultados[,c("SYMBOL","TRANSCRIPT", "ACCESSION","REFSEQ_ID",
                                    "ENTREZ_GENE_ID", "PROBE_ID","CHROMOSOME",
                                    "PROBE_CHR_ORIENTATION",
                                    "logFC","P.Value","adj.P.Val","B")]

volcanoplot(PAREfit2, coef="PAREgrupos2", highlight=10, names=PAREfit2$genes$SYMBOL) #mudar

PAREtop100_FOCA <- topTable(PAREfit2, coef="PAREgrupos2", adjust="bonferroni", 
                            number=100) #mudar
PAREtop100_FOCA <- PAREtop100_FOCA[,c("SYMBOL","TRANSCRIPT", "ACCESSION","REFSEQ_ID",
                                      "ENTREZ_GENE_ID", "PROBE_ID","CHROMOSOME",
                                      "PROBE_CHR_ORIENTATION",
                                      "logFC","P.Value","adj.P.Val","B")]

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:

write.table(PAREtop100_FOCA,file="PAREtop100_BaixoDesvioCBCL_coca.txt", sep="\t", 
            row.names=F,quote=FALSE) #mudar
write.table(PAREresultados,file="PAREresultados_BaixoDesvioCBCL_coca.txt", sep="\t", 
            row.names=F,quote=FALSE) #mudar
################################################################
#FIM DOs Pareados baixo desvio, curso normal n˜ão relacionado PSY?
###############################################################

################################################################
#INICIO DOs Pareado alto desvio
###############################################################
# identificar as duplas e o grupo de comparação

PAREgrupos <- factor(phenoData(eset_PAREbaixodesvio)$GROUPS, levels=c("1", "2")) #mudar
ID_INPD <- factor(phenoData(eset_PAREbaixodesvio)$Sample_ID) #mudar

#GedBlood
# PAREdesign <- model.matrix(~ID_PEP+PAREgrupos+neutro+NK+monogb+igg+igm+tc+tc_act)
PAREdesign <- model.matrix(~ID_INPD+PAREgrupos)
PAREdesign <- model.matrix(~PAREgrupos+ID_INPD)

PAREfit <- lmFit(eset_PAREbaixodesvio, PAREdesign) #mudar

PAREfit2 <- eBayes(PAREfit)
# PAREresults <- decideTests(PAREfit2)
# PAREresults2 <- decideTests(PAREfit2, adjust="bonferroni")

PAREresultados <- topTable(PAREfit2,coef="PAREgrupos2", adjust="bonferroni", 
                           number=Inf) #mudar

topTable(PAREfit2,coef="PAREgrupos2", adjust="bonferroni", 
         number=Inf) #mudar

PAREresultados <- PAREresultados[,c("SYMBOL","TRANSCRIPT", "ACCESSION","REFSEQ_ID",
                                    "ENTREZ_GENE_ID", "PROBE_ID","CHROMOSOME",
                                    "PROBE_CHR_ORIENTATION",
                                    "logFC","P.Value","adj.P.Val","B")]

volcanoplot(PAREfit2, coef="PAREgrupos2", highlight=10, names=PAREfit2$genes$SYMBOL) #mudar

PAREtop100_FOCA <- topTable(PAREfit2, coef="PAREgrupos2", adjust="bonferroni", 
                            number=100) #mudar
PAREtop100_FOCA <- PAREtop100_FOCA[,c("SYMBOL","TRANSCRIPT", "ACCESSION","REFSEQ_ID",
                                      "ENTREZ_GENE_ID", "PROBE_ID","CHROMOSOME",
                                      "PROBE_CHR_ORIENTATION",
                                      "logFC","P.Value","adj.P.Val","B")]

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:

write.table(PAREtop100_FOCA,file="PAREtop100_BaixoDesvioCBCL_coca.txt", sep="\t", 
            row.names=F,quote=FALSE) #mudar
write.table(PAREresultados,file="PAREresultados_BaixoDesvioCBCL_coca.txt", sep="\t", 
            row.names=F,quote=FALSE) #mudar

###############################################################
#FIM DOs Pareado alto desvio
###############################################################

###############################################################
#INICIO DO Any disorder
###############################################################
# Variáveis que vou utilizar e que estão no meu pDat
row.names(juntado) <- juntado$sampleID
teta <- eset_ANYDIS_1w@phenoData@data
juntado2 <- juntado[teta$sampleID,]

grupos <- factor(juntado2$dcFUPany, levels=c("0", "2")) #mudar
sex <- factor(phenoData(eset_ANYDIS_1w)$SEX, levels=c("MALE", "FEMALE")) #mudar
# PC1 <- phenoData(eset_caco_0w)$PC1_all_pre_ComBat

####################################
# quais serão as minhas covariaveis? 
####################################
#GedBlood
# Este foi utilizado o artigo 4
design <- model.matrix(~grupos+sex) 

# design <- model.matrix(~0+grupos)

#Fit
fit <- lmFit(eset_ANYDIS_1w, design) #mudar
fit2 <- contrasts.fit(fit, coef="grupos2") #mudar
fit3 <- eBayes(fit2)
topTable(fit3, coef="grupos2", adjust="bonferroni") #mudar
results <- decideTests(fit2)
results_bonferroni <- decideTests(fit2, adjust="bonferroni")

# Visualizar resultados com FDR e com bonferroni
#FDR
vennDiagram(results)
vennDiagram(results, include="up")
vennDiagram(results, include="down")
#Bonferroni
vennDiagram(results_bonferroni)
vennDiagram(results_bonferroni, include="up")
vennDiagram(results_bonferroni, include="down")

#hora de organizar os meus resultados
resultados_CACO <- topTable(fit3, coef="grupos2", adjust.method="bonferroni", 
                            number=Inf) #mudar

#volcano plot com o resultado de eBayes
volcanoplot(fit3, coef="grupos2", highlight=20, names=fit3$genes$SYMBOL) #mudar

#Só os top 100
#top1000_CACO <- topTable(fit3, coef="gruposCASE", adjust.method="bonferroni", number=1000)
top100_CACO <- topTable(fit3, coef="grupos2", adjust.method="bonferroni", number=100) #mudar

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:
write.table(resultados_CACO, file="Resultados_ANYDISORDER.txt",sep="\t", 
            row.names=F,quote=FALSE)
write.table(top100_CACO,file="top100_ANYDISORDER.txt", sep="\t", row.names=F,
            quote=FALSE)

###############################################################
#Fim  DO ANY DISORDER
###############################################################

###############################################################
#INICIO DO Adrugs
###############################################################
# Variáveis que vou utilizar e que estão no meu pDat
row.names(juntado2) <- juntado2$sampleID
teta <- eset_drug_1w@phenoData@data
juntado3 <- juntado2[teta$sampleID,]

grupos <- factor(juntado3$RiskFUPASUseever, levels=c("0", "1")) #mudar
sex <- factor(phenoData(eset_drug_1w)$SEX, levels=c("MALE", "FEMALE")) #mudar
# PC1 <- phenoData(eset_caco_0w)$PC1_all_pre_ComBat

####################################
# quais serão as minhas covariaveis? 
####################################
#GedBlood
# Este foi utilizado o artigo 4
design <- model.matrix(~grupos+sex) 

# design <- model.matrix(~0+grupos)

#Fit
fit <- lmFit(eset_drug_1w, design) #mudar
fit2 <- contrasts.fit(fit, coef="grupos1") #mudar
fit3 <- eBayes(fit2)
topTable(fit3, coef="grupos1", adjust="bonferroni") #mudar
results <- decideTests(fit2)
results_bonferroni <- decideTests(fit2, adjust="bonferroni")

# Visualizar resultados com FDR e com bonferroni
#FDR
vennDiagram(results)
vennDiagram(results, include="up")
vennDiagram(results, include="down")
#Bonferroni
vennDiagram(results_bonferroni)
vennDiagram(results_bonferroni, include="up")
vennDiagram(results_bonferroni, include="down")

#hora de organizar os meus resultados
resultados_CACO <- topTable(fit3, coef="grupos1", adjust.method="bonferroni", 
                            number=Inf) #mudar

#volcano plot com o resultado de eBayes
volcanoplot(fit3, coef="grupos1", highlight=20, names=fit3$genes$SYMBOL) #mudar

#Só os top 100
#top1000_CACO <- topTable(fit3, coef="gruposCASE", adjust.method="bonferroni", number=1000)
top100_CACO <- topTable(fit3, coef="grupos1", adjust.method="bonferroni", number=100) #mudar

#############
# ATENÇÃO!!! Checar o nome que vai salvar
#############
#Salvando:
write.table(resultados_CACO, file="Resultados_Drug.txt",sep="\t", 
            row.names=F,quote=FALSE)
write.table(top100_CACO,file="top100_Drug.txt", sep="\t", row.names=F,
            quote=FALSE)
