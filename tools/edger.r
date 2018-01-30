library(optparse)
library(edgeR)
library(methods)

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
dds <- DGEList(counts=count_mtx, group=factor(condition))
design_mtx <- model.matrix(~0 + dds$samples$group)
colnames(design_mtx) <- levels(dds$samples$group)

# ## run EdgeR
# cat('... running EdgeR\n')
dds <- calcNormFactors(dds)
# dds <- estimateDisp(dds, design_mtx)
# dds <- estimateGLMCommonDisp(dds, design_mtx)
dds <- estimateGLMTrendedDisp(dds, design_mtx)
dds <- estimateGLMTagwiseDisp(dds, design_mtx)
fit <- glmFit(dds, design_mtx)
lrt <- glmLRT(fit)
res <- topTags(lrt, n=dim(count_mtx)[1], adjust.method="BH", sort.by="PValue")

# dds <- calcNormFactors(dds)
# dds <- estimateCommonDisp(dds)
# dds <- estimateTagwiseDisp(dds)
# de <- exactTest(dds, pair=c("untreated", "treated"))
# res <- topTags(de, n=dim(count_mtx)[1], adjust.method="BH", sort.by="PValue")

## write result
cat('... writing output\n')
write.table(res, file=opt$file_output, quote=FALSE, sep='\t')