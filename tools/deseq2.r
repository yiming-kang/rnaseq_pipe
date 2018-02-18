library(optparse)
library(xlsx)
# library(DESeq2)

parse_arguments <- function(){
	option_list <- list(
		make_option(c('-c', '--count_mtx'), 
					help='Read count matrix (genes x samples).'),
		make_option(c('-d', '--design_table'), 
					help='Design table containing sample grouping indicator.'),
		make_option(c('-q', '--qa_table'), 
					help='QC table containing audit stauts.'),
		make_option(c('-o', '--output_dir'), 
					help='Table of DE genes ranked by ajusted p-value.'))
	opt <- parse_args(OptionParser(option_list=option_list))
	return(opt)
}


parse_metadata <- function(design_filepath, qa_filepath) {
	## Parse design and QA tables
	design <- read.xlsx(design_filepath, sheetIndex=1, check.names=FALSE)
	print(colnames(design))
	q()
}


## main
parsed_opt <- parse_arguments()
parse_metadata(parsed_opt$design_table, parsed_opt$qa_table)

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