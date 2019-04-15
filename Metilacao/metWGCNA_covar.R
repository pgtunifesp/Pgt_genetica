##################################
#FASE 1- Preparação do dado
###################################

## Chamar as bibliotecas

library(WGCNA)
library(methylumi)
library(TCGAMethylation450k)
library(wateRmelon)
library(minfi)
require(plyr)

# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "~/Documents/PEP/Methylation/Data/";
setwd(workingDir)

#vou chamar os dados clínicos primeiro, diferente do tutorial para selecionar casos e controles
traitData <- read.csv("~/Downloads/BR_FEP_Meth_samplesheet_final290116.csv", header=T)
dim(traitData)
names(traitData)

##Minu: função para multisub
mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}

# vars e covars que ficam
vars <- c("Sample_Name","PHENOTYPE","SEX","AGE","SMOKE","SMOKED","COR","FUMANTE","EXPRE","barcodes")
# vars <- c("Sample_Name","PHENOTYPE","barcodes")
#vars <- c("sampleID","PHENOTYPE","SEX","AGE")
#vars <- c("sampleID", "PHENOTYPE","SEX","AGE","SMOKE")


# remove columns that hold information we do not need.
allTraits = traitData[, vars];
fock <- allTraits$PHENOTYPE == "FEP" | allTraits$PHENOTYPE == "WFOLLOW"
FOCA <- allTraits[fock,]

FOCA$PHENOTYPE <- mgsub(c("FEP","WFOLLOW"),c(1,2),FOCA$PHENOTYPE)
FOCA$SEX <- mgsub(c("FEMALE","MALE"),c(2,1),FOCA$SEX)

head(FOCA)
dim(FOCA)
names(FOCA)

#Chamar os resultados do pipeline foca2_corrigido
gc()
beta <- read.csv("./Betanias_dez17.csv", header = T)
dim(beta)
beta[1:5,1:5]

fock2 <- names(beta) %in% FOCA$barcodes
beta_FOCA <- beta[,fock2]
# write.table(beta_FOCA,"betasFOCA.txt", quote = F,row.names = T,col.names = T)
dim(beta_FOCA)
#
#Precisa tranformar em character para reconhecer como coluna depois...
#allTraits$Sample_Name <- as.character(allTraits$Sample_Name)
#allTraits$barcodes <- as.character(allTraits$barcodes)
#IDs <- allTraits$barcodes
#beta2 <- beta[,IDs]
#escreva o meu csv, mas sou malandro, chame de tabela mas separ com virgula e tira as aspas pq elas dão treta
#write.table(mutreta4, file="mutretagem.csv", quote=F, row.names=T, col.names=T, sep=",")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the FEP BR data set
datExpr0 = as.data.frame(t(beta_FOCA));

head(datExpr0)[1:5,1:5]

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
# checar quem são eles
abline(h = 250, col = "red");
# Determine cluster under the line
#minu:pilantrei
clust = cutreeStatic(sampleTree, cutHeight = 250, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Minu: Tirei a chamada do banco clinico daqui...

#Incluir o X para ficar igual ao datExp
#Xis <- rep("X",length(traitData$Study_ID))
#traitData2 <- cbind(Xis, traitData)
#traitData2$ID <-paste(Xis, traitData2$sampleID, sep="")
#dim(traitData2)


# Form a data frame analogous to expression data that will hold the clinical traits.
FEPSamples = rownames(datExpr);
traitRows = match(FEPSamples, FOCA$barcodes);
#Se funcionou bateu o ID de todos
traitRows
datTraits = FOCA[traitRows, -1];
rownames(datTraits) = FOCA[traitRows, 1];
head(datTraits)
collectGarbage();

#Elas precisam ser numericas
datTraits$SEX <- as.numeric(as.numeric(datTraits$SEX))
datTraits$PHENOTYPE <- as.numeric(as.character(datTraits$PHENOTYPE))
datTraits$AGE <- as.numeric(as.character(datTraits$AGE))
datTraits$SMOKE <- as.numeric(as.character(datTraits$SMOKE))
datTraits$SMOKED <- as.numeric(as.character(datTraits$SMOKED))
datTraits$FUMANTE <- as.numeric(as.character(datTraits$FUMANTE))
datTraits$COR <- mgsub(c("branca","amarela", "mulato","negra","outra"),c(1,2,3,4,5),datTraits$COR)
datTraits$COR <- as.numeric(as.character(datTraits$COR))
# datTraits$NK <- as.numeric(as.character(datTraits$NK))
# datTraits$TC <- as.numeric(as.character(datTraits$TC))
# datTraits$TC_ACT <- as.numeric(as.character(datTraits$TC_ACT))
# datTraits$IGG <- as.numeric(as.character(datTraits$IGG))
# datTraits$IGM <- as.numeric(as.character(datTraits$IGM))
# datTraits$DC_AC <- as.numeric(as.character(datTraits$DC_AC))
# 
fica <- c(1:8)
# 
datTraits <- datTraits[,fica]

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file = "FEPBR-01-dataInput.RData")


###################################
#FASE 2!!!!!!!
#Network construction and module detection
#Opção a: Automatic, one-step network construction and module detection
###################################
#=======================Met_vars <- subset(SD, SD > 0.01)
==============================================================
#
#  Code chunk 1
#
#=====================================================================================
workingDir = "~/Documents/PEP/Methylation/Data/";
setwd(workingDir)
library(WGCNA)
disableWGCNAThreads()
# enableWGCNAThreads(nThreads=2)
# Load the data saved in the first part
lnames = load(file = "FEPBR-01-dataInput.RData");

#The variable lnames contains the names of loaded variables.
lnames

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
# sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# # Plot the results:
# sizeGrWindow(9, 5)
# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
softpower = 6


adjacency = adjacency(datExpr,power = softpower)

# Precisava olhar o grafico e mudar de acordo com a reta, no meu caso mudei para 8.

net = blockwiseModules(datExpr, power = 6, networkType = "unsigned",
                                TOMType = "unsigned",
                                saveTOMs = TRUE,
                                saveTOMFileBase = "FEPMouseTOM", 
                                verbose = 3, maxBlockSize = 20000,minModuleSize = 1000)

net_unsigned = blockwiseModules(datExpr, power = 9, networkType = "unsigned",
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "FEPMouseTOM", 
                       verbose = 3)
# Quantos blocos eu tenho, e quantos genes em cada bloco, o bloco ) eh para aqueles
# genes solitários e sem ninguem.
table(net_unsigned$colors)
table(net$colors)

net = net_supersigned

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.1,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "FEPBR-02-networkConstruction-auto.RData")


###################################
#FASE3!!
#Relating modules to external clinical traits and identifying important genes
###################################

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================
heatmap

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

# Define variable pheno containing the pheno column of datTrait
#mudei mas não sei o que eh
pheno = as.data.frame(datTraits$PHENOTYPE);
names(pheno) = "Phenotype"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, pheno, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(pheno), sep="");
names(GSPvalue) = paste("p.GS.", names(pheno), sep="");


#=====================================================================================
#
#  Code chunk 5
###IMPORTANTEEEE, PRECISA MUDAR A SUA COR!!
#=====================================================================================

module = "pink"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body pheno",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================

str(eset_final)

nuID <-eset_final@featureData@data

checagem <- cbind(names(datExpr),nuID)
head(checagem)
tail(checagem)


colnames(datExpr) <- checagem$ENTREZ_GENE_ID

names(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================

# Isso tem que ser feito para cada modulo significante!

red <- names(datExpr)[moduleColors=="red"]
tan <- names(datExpr)[moduleColors=="tan"]
turquoise <- names(datExpr)[moduleColors=="turquoise"]
yellow <- names(datExpr)[moduleColors=="yellow"]
purple <- names(datExpr)[moduleColors=="purple"]
#module6 <- names(datExpr)[moduleColors=="green"]
#module7 <- names(datExpr)[moduleColors=="green"]

intModules <- c("purple","red","tan","turquoise","yellow")


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================

probes = names(datExpr)

#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================
# Nesse passo precisa mudar 
# Create the starting data frame
geneInfo0 = data.frame(RefSeqID = probes,
                       geneSymbol = checagem$SYMBOL,
                       LocusLinkID = checagem$ENTREZ_GENE_ID,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for pheno
modOrder = order(-abs(cor(MEs, pheno, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Phenotype));
geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================

write.csv(geneInfo, file = "geneInfo.csv")

#############################
#FASE 4!!!
#Interfacing network analysis with other data such as functional annotation and gene ontology
#############################
#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================
#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Read in the probe annotation
#annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
probes = names(datExpr)
probes2annot = match(probes, checagem$ENTREZ_GENE_ID)
# Get the corresponding Locuis Link IDs
allLLIDs = checagem$ENTREZ_GENE_ID[probes2annot];
# $ Choose interesting modules
#MINU:
intModules
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("ENTREZ_GENE_IDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("ENTREZ_GENE_IDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "human", nBestP = 10);

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

tab = GOenr$bestPTerms[[4]]$enrichment

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

names(tab)

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================

write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)

#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================

keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab

########################
#FASE 5!!
#Network visualization using WGCNA functions
########################

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================
#Minu: deu erro...

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 9);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)

TOMplot(plotTOM, geneTree, moduleColors)

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

nSelect = 4000
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate pheno from the clinical traits
pheno = as.data.frame(datTraits$PHENOTYPE);
names(pheno) = "pheno"
# Add the pheno to existing module eigengenes
MET = orderMEs(cbind(MEs, pheno))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


################
#FASE 6!!!!
#Export of networks to external software
################

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
#annot = read.csv(file = "GeneAnnotation.csv");
#Teste!!!
intModules 
# 
for (module in intModules)
{ 
  # Select module probes
  inModule = (moduleColors==module)
  # 
  modProbes = probes[inModule];
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into an edge list file VisANT can read
  vis = exportNetworkToVisANT(modTOM,
                              file = paste("VisANTInput-", module, ".txt", sep=""),
                              weighted = TRUE,
                              threshold = 0,
                              probeToGene = data.frame(checagem$ENTREZ_GENE_ID, checagem$SYMBOL) )
  
  
}

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

nTop = 30;
IMConn = softConnectivity(datExpr);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(TOM[top, top],
                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(checagem$REFSEQ_ID, checagem$SYMBOL) )

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
#annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
modules = intModules
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = checagem$SYMBOL[match(modProbes, checagem$ENTREZ_GENE_ID)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);
