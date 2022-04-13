library("dplyr")
library("DESeq2")

counts = snakemake@input[['counts']]

metadata <- snakemake@params[['samples']]

sampleID <- snakemake@params[['sample_id']]

Type <- snakemake@params[['linear_model']]

outCounts = snakemake@output[["counts"]]

outLogCounts = snakemake@output[["logcounts"]]

# Read in metadata table and order according to sampleID
md <- read.delim(file=metadata, sep = "\t", stringsAsFactors = FALSE)
md <- md[order(md[sampleID]),]

# Read in counts table
subdata <- read.table(counts, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
subdata <- subdata[,order(colnames(subdata))]

# Keep only the PP_IDs of the types we have chosen in the metadata table above
rownames(md) <- md[[sampleID]]
md[[sampleID]] <- NULL
keep <- colnames(subdata)[colnames(subdata) %in% rownames(md)]
subdata <- subdata[, keep]
dim(subdata)

# Check
stopifnot(rownames(md)==colnames(subdata))

# Obtain the number of genes that meet padj<0.01 for reference line in histogram
dds <- DESeqDataSetFromMatrix(countData=subdata,
                              colData=md,
                              design= as.formula(paste('~',Type)))

dds <- estimateSizeFactors(dds)

# Remove uninformative columns
dds <- dds[ rowSums(counts(dds)) >= 1, ]

# Normalization and pre-processing
dds <- DESeq(dds)

dds <- dds[ rowSums(counts(dds)) > 1, ]

write.table(counts(dds, normalized = TRUE), outCounts, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(log2(counts(dds, normalized = TRUE)), outLogCounts, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
