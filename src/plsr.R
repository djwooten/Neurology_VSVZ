library(CMA)

cross.validation.stats <- function(cv) {

    l = length(cv)
    df = data.frame(row.names=1:l)
    df['missclassifications']=NA
    df['sensitivity']=NA
    df['specificity']=NA
    for (ind in 1:l) {
        x = cv[[ind]]

        dims <- sort(unique(c(x@y, x@yhat))) 
        ftab <- matrix(0, nrow=length(dims), ncol=length(dims))
        for(i in seq(along=dims)){
            for(j in seq(along=dims)){
                ftab[i,j] <- length(intersect(which(x@y==dims[i]), which(x@yhat == dims[j])))
            }
        } 
        ftab <- as.table(ftab)

        df$missclassifications[ind] = (sum(ftab)-sum(diag(ftab)))/sum(ftab)
        if(x@mode == "binary"){
            df$sensitivity[ind] = ftab[2,2]/sum(ftab[2,])
            df$specificity[ind] = ftab[1,1]/sum(ftab[1,])
        }
    }

    return(df)
    
}

plot.cv <- function(df,fname="cv.pdf",l=-3) {
    pdf(fname,height=6,width=12)
    par(mfrow=c(1,3), cex.lab=1.5)
    with(df, hist(missclassifications,xlim=c(0,1),xlab="Missclassification Rate",main="",breaks=seq(0,1,0.025)))
    with(df, hist(sensitivity,xlim=c(0,1),xlab="Sensitivity",main="",breaks=seq(0,1,0.025),ylab=""))
    with(df, hist(specificity,xlim=c(0,1),xlab="Specificity",main="",breaks=seq(0,1,0.025),ylab=""))
    mtext("PLS-LR cross validation accuracy statistics (n=100)", side = 3, line = l, outer = TRUE, cex=1.5)
    dev.off()
}

data = read.csv("../expression.csv",row.names=1,na.strings=c("null","na","NA","NULL","nan","NaN","NAN"))

patient.cols <- colnames(data)
for (i in 1:length(patient.cols)) { patient.cols[i] <- strsplit(as.character(patient.cols[i]),'\\.')[[1]][3] }

clusters <- read.csv("../data/patient_status.csv",row.names=1)
clusters.rows <- row.names(clusters)
for (i in 1:length(clusters.rows)) { clusters.rows[i] <- strsplit(as.character(clusters.rows[i]),'-')[[1]][3] }

svz.has <- patient.cols %in% clusters.rows
svz.order <- match(patient.cols[svz.has],clusters.rows)
t
data.svz <- data[,svz.has]
colnames(data.svz) <- patient.cols[svz.has]

 
svz <- clusters[svz.order,'SVZ']


ratio <- 2/3
set.seed(111)
data.svz.t <- t(data.svz)

learning_sets <- GenerateLearningsets(y=row.names(data.svz.t),method="MCCV",ntrain=floor(ratio*length(row.names(data.svz.t))),niter=100)

cma.result.cl <- classification(data.svz.t, svz, learningsets=learning_sets, classifier=pls_lrCMA)

pdf("plsr_prediction.pdf")
plot(cma.result.cl[[1]])
dev.off()

df <- cross.validation.stats(cma.result.cl)
plot.cv(df)


