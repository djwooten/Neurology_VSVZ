library(ConsensusClusterPlus)
library(gplots)
library(devtools)

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

inputfileDir <- '/data/dwooten/School/qlab/Projects/GBM/ht_hg_u133a'
workingDir <- getwd()

# Read data
print("Loading Data")
setwd(inputfileDir)
data <- read.csv('expression.csv',row.names=1)
setwd(workingDir)

patient.cols <- colnames(data)
for (i in 1:length(patient.cols)) { patient.cols[i] <- strsplit(as.character(patient.cols[i]),'\\.')[[1]][3] }

clusters <- read.csv("../data/patient_status.csv",row.names=1)
clusters.rows <- row.names(clusters)
for (i in 1:length(clusters.rows)) { clusters.rows[i] <- strsplit(as.character(clusters.rows[i]),'-')[[1]][3] }

subtype.clusters <- read.csv("../data/patient_status_subtype.csv",row.names=1)
subtype.clusters.rows <- row.names(subtype.clusters)
for (i in 1:length(subtype.clusters.rows)) { subtype.clusters.rows[i] <- strsplit(as.character(subtype.clusters.rows[i]),'-')[[1]][3] }

# Median center by gene
print("Median center")
data <- sweep(data,1,apply(data,1,median,na.rm=T))

# Function for plotting heatmap + colorbars
plot.heatmap <- function(r, k) {
    if (k>4) return
    mydist <- function(c) {0}
    myclust <- function(c) {r$consensusTree}

    #vsvz.colors <- clusters[match(patient.cols[r$consensusTree$order],clusters.rows),"SVZ"]
    vsvz.colors <- clusters[match(patient.cols,clusters.rows),"SVZ"]
    vsvz.colors[vsvz.colors==0] <- 'black'
    vsvz.colors[vsvz.colors==1] <- 'yellow'
    vsvz.colors[is.na(vsvz.colors)] <- 'gray60'

    #subtype.colors <- clusters[match(patient.cols[r$consensusTree$order],clusters.rows),"Subtype"]
    subtype.colors <- as.character(subtype.clusters[match(patient.cols,subtype.clusters.rows),"Subtype"])
    subtype.colors[subtype.colors=="Neural"] <- 'coral'
    subtype.colors[subtype.colors=="Proneural"] <- 'orangered4'
    subtype.colors[subtype.colors=="Mesenchymal"] <- 'darkgreen'
    subtype.colors[subtype.colors=="G-CIMP"] <- 'turquoise1'
    subtype.colors[subtype.colors=="Classical"] <- 'green'
    subtype.colors[subtype.colors=="N/A"] <- 'gray60'
    subtype.colors[is.na(subtype.colors)] <- 'gray60'

    #cc.colors <- as.numeric(r$consensusClass)[r$consensusTree$order]
    cc.colors <- as.numeric(r$consensusClass)
    cc.colors[cc.colors==1] <- 'darkorchid'
    cc.colors[cc.colors==2] <- 'hotpink'
    cc.colors[cc.colors==3] <- 'gold4'
    cc.colors[cc.colors==4] <- 'midnightblue'
    cc.colors[is.na(cc.colors)] <- 'gray60'

    legend.labels <- c("VSVZ-",
                       "VSVZ+",
                       "",########################
                       "CC 1",
                       "CC 2", 
                       "CC 3", 
                       "CC 4", 
                       "", ########################
                       "Neural",
                       "Proneural",
                       "Mesenchymal",
                       "G-CIMP",
                       "Classical", 
                       "", ########################
                       "N/A")



    legend.colors <- c("black",
                       "yellow",
                       "white",########################
                       "darkorchid",
                       "hotpink",
                       "gold4",
                       "midnightblue",
                       "white",########################
                       "coral",
                       "orangered4",
                       "darkgreen",
                       "turquoise1",
                       "green",
                       "white",########################
                       "gray60")

    if (k == 2) {
        legend.labels <- legend.labels[-c(6,7)]
        legend.colors <- legend.colors[-c(6,7)]
    }
    if (k == 3) {
        legend.labels <- legend.labels[-c(7)]
        legend.colors <- legend.colors[-c(7)]
    }

    clab <- cbind(subtype.colors, cc.colors, vsvz.colors)
    colnames(clab) <- c("Subtype", "Consensus Cluster", "VSVZ Status")

    myPal = function(n=10){
        #returns n colors
        seq = rev(seq(0,255,by=255/(n)))
        palRGB = cbind(seq,seq,255)
        rgb(palRGB,maxColorValue=255)
    }

    pdf(file=paste(k,'_clusters.pdf',sep=""))
    
    main_title="Consensus Clusters"
    par(cex.main=1)

    heatmap.3(r$consensusMatrix, 
     hclustfun=myclust, 
     distfun=mydist, 
     na.rm = TRUE, 
     scale="none", 
     dendrogram="column", 
     margins=c(6,12),
     Rowv=as.dendrogram(r$consensusTree), Colv=as.dendrogram(r$consensusTree),
     ColSideColors=clab, 
     symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", 
     main=main_title, 
     labCol=FALSE, labRow=FALSE, 
     #col=rev(heat.colors(75)),
     col=myPal(),
     ColSideColorsSize=3, RowSideColorsSize=0, 
     KeyValueName="Consensus Index")
     legend("topright", legend=legend.labels, fill=legend.colors, border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
    dev.off()
}




clusterAlg = 'kmdist'
distance = 'pearson
outputfileDir <- paste(clusterAlg,distance,sep="_")

if (!dir.exists(outputfileDir)) dir.create(outputfileDir)


# Consensus clustering
print("Consensus clustering")
results <- ConsensusClusterPlus(as.matrix(data)
                                ,maxK        =   10
                                ,reps        =   1000
                                ,pItem       =   0.8
                                ,pFeature    =   0.8
                                ,title       =   outputfileDir
                                ,clusterAlg  =   clusterAlg
                                ,distance    =   distance
                                ,seed        =   1
                                ,plot        =   "png"
                                )

setwd(outputfileDir)
for (k in c(2:10)) {
    lbls <- results[[k]]$consensusClass
    treeOrder <- results[[k]]$consensusTree$order
    df <- data.frame(row.names=names(lbls)[treeOrder])
    df$class <- as.numeric(lbls)[treeOrder]
    write.table(df, file=paste(k,'_clusters.txt',sep=""), col.names=FALSE, quote=FALSE)
    plot.heatmap(results[[k]],k)
}
saveRDS(results,"results.rds")
setwd(workingDir)
}
