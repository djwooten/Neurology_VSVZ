library(WGCNA)

# Set options
options(stringsAsFactors = FALSE)
allowWGCNAThreads()


###############################################
# STEP 1: LOAD DATA
###############################################
# Load data
data_exp = read.csv("../gbm_expression.csv",row.names=1)
data = as.data.frame(t(data_exp))

nGenes = ncol(data)
nSamples = nrow(data)

###############################################
# STEP 2: NETWORK CONSTRUCTION
###############################################

# Choose soft-threshold to maximize scale-free network topology
powers = c(c(1:10), seq(from=12, to=20, by=2))
sft = pickSoftThreshold(data, powerVector=powers, verbose=5, networkType="signed")
sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1 = 0.9

##########################################
# STEP 3: MANUAL INTERVENTION TO IDENTIFY POWER
##########################################
pdf("powers_signed.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
    main=paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers, cex=cex1, col="red")
abline(h=0.9, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main=paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

softPower = 12 # Manually determined from powers_signed.pdf


##########################################
# STEP 4: MODULE CONSTRUCTION
##########################################

# Compute adjacency, TOM
adjacency = adjacency(data, power=softPower, type="signed")
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Make clusters
geneTree = hclust(as.dist(dissTOM), method="average")
minModuleSize = 100
dynamicMods = cutreeDynamic(dendro=geneTree, distM = dissTOM, deepSplit=2, pamRespectsDendro=FALSE, minClusterSize=minModuleSize)
dynamicColors = labels2colors(dynamicMods)

pdf("pwr12/signed_WGCNA_CCLE.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
    dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05,
    main="Gene dendrogram and module colors")
dev.off()

table(dynamicColors)
save(geneTree, dynamicColors, file="pwr12/trees_and_colors.RData")

for (which.module in unique(dynamicColors)) {
    module.genes <- names(data[,dynamicColors==which.module])
    write.table(module.genes, file=paste("pwr12/genes_",which.module,".txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
}

MEs =  moduleEigengenes(data, colors = dynamicColors)$eigengenes
rownames(MEs) <- rownames(data)
write.csv(MEs, file="./pwr12/module_eigengenes.csv")
kMEs = signedKME(data, MEs)
write.csv(kMEs, file="./pwr12/kMEs.csv")
