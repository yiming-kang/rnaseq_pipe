library(optparse)
library(DESeq2)

option_list <- list(
	make_option(c('-i', '--file_count_mtx'), help='Read count matrix (genes x samples)'),
	make_option(c('-o', '--file_output'), help='Table of DE genes ranked by ajusted p-value'))
opt <- parse_args(OptionParser(option_list=option_list))

## parse count matrix in data.frame format
cat('... preparing data\n')
count_mtx <- read.csv(opt$file_count_mtx, check.names=FALSE, row.names=1)
samples <- colnames(count_mtx)
genes <- rownames(count_mtx)

## prepare count and design matrices for DESeq2
cat('... preparing design matrix\n')
condition <- c(rep('untreated',4), rep('treated',4))
coldata <- data.frame(condition, row.names=samples)
dds <- DESeqDataSetFromMatrix(countData=count_mtx, colData=coldata, design=~condition)

## run DESeq2
cat('... running DESeq2\n')
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]

## write result
cat('... writing output\n')
write.table(res, file=opt$file_output, quote=FALSE, sep='\t')