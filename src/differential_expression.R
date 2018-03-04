library(nortest)

set.seed(120424)


#data = read.csv("../expression.csv",row.names=1,na.strings=c("null","na","NA","NULL"))

# Remove genes with NA (or null)
data.cc <- data[complete.cases(data), ]

# Strip out a list of the 4-digit patient IDs for every sample
patient.cols <- colnames(data)
for (i in 1:length(patient.cols)) { patient.cols[i] <- strsplit(as.character(patient.cols[i]),'\\.')[[1]][3] }

# Read SVZ status
clusters <- read.csv("../data/patient_status.csv",row.names=1)
clusters.rows <- row.names(clusters)

# Reduce each patient name to the 4-digit patient ID
for (i in 1:length(clusters.rows)) { clusters.rows[i] <- strsplit(as.character(clusters.rows[i]),'-')[[1]][3] }

# Find set of patients that have 0 or 1 SVZ status
svz.has <- patient.cols %in% clusters.rows

# Match: for each value in patient.cols[svz.has], find its index in clusters.rows
svz.order <- match(patient.cols[svz.has],clusters.rows)

# Subset to just that dataset
data.svz <- data[,svz.has]
colnames(data.svz) <- patient.cols[svz.has]

# Read the svz status of those
svz <- clusters[svz.order,'SVZ']

# Compute statistical test between patients with svz==0 and svz==1
resut <- apply(data.svz, 1, function(x) wilcox.test(as.numeric(x[svz==0]),as.numeric(x[svz==1]),paired=FALSE,na.rm=TRUE)$p.value)
#result <- apply(data.svz, 1, function(x) t.test(as.numeric(x[svz==0]),as.numeric(x[svz==1]),paired=FALSE,na.rm=TRUE)$p.value)

# Do FDR correction
result.bh <- p.adjust(result,method="BH")

# Statistical test results put in dataframe to be saved
df <- data.frame(row.names=names(result))
df$p <- result
df$bh <- result.bh
write.csv(df,"wilcox_pvals_eigengenes.csv")


# Test normality
normality <- data.frame(row.names=names(result))
normality$svz0 <- apply(data.svz,1,function(x) ad.test(as.numeric(x[svz==0]))$p.value)
normality$svz1 <- apply(data.svz,1,function(x) ad.test(as.numeric(x[svz==1]))$p.value)
write.csv(normality,"normality.csv")


