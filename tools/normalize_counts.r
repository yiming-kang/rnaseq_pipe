suppressMessages(library(optparse))
suppressMessages(library(DESeq2))


parse_arguments <- function() {
	option_list <- list(
		make_option(c('-i', '--file_count_mtx'), 
					help='Read count matrix (genes x samples).'),
		make_option(c('-o', '--output'), 
					help='Filepath or directory path for normalized count file(s).'))
	args <- parse_args(OptionParser(option_list=option_list))
	return(args)
}


read_count_matrix <- function(filepath) {
	## read in count matrix
	count_mtx <- read.csv(filepath, check.names=FALSE, row.names=1)
	count_mtx <- count_mtx[order(row.names(count_mtx)),]
	return(count_mtx)
}

normalize_count <- function(count_mtx, samples) {
	## normalize with size factor
	n_samples <- length(samples)
	conditions <- c(rep('dummy0',floor(n_samples/2)), 
					rep('dummy1',n_samples-floor(n_samples/2)))
	coldata <- data.frame(conditions, row.names=samples)
	dds <- DESeqDataSetFromMatrix(countData=count_mtx, colData=coldata, 
									design=~conditions)
	dds <- estimateSizeFactors(dds)
	norm_count_mtx <- counts(dds, normalized=TRUE)
	return(norm_count_mtx)
}


write_count_matrix_file <- function(filepath_output, count_mtx) {
	## write count matrix file
	write.csv(count_mtx, file=filepath_output, quote=FALSE)
}


write_indiv_count_file <- function(dirpath_output, count_mtx, samples, gene_lengths=NULL) {
	## write output files
	for (sample in samples) {
		filepath <- paste(dirpath_output, sample, sep='/')
		count_out <- count_mtx[,c(sample)]
		write.table(count_out, file=paste0(filepath,'.count'), 
					quote=FALSE, row.names=FALSE, col.names=FALSE)
		# ## caluclate FPKMs if gene lengths are given
		# if (!is.null(gene_lengths)) {
		# 	fpkm_out <- count_out / (genes_len*sum(count_out)) * 10^9
		# 	write.table(fpkm_out, file=paste0(filepath,'.fpkm'), quote=FALSE, row.names=FALSE, col.names=FALSE)
		# }
	}
}


# get_gene_lengths <- function(genes) {
# 	## load gene length dictionary
# 	annot <- read.table('H99/tmp', sep='\t', check.names=FALSE, header=TRUE, comment.char='@')
# 	gene_len_dict <- list()
# 	for (gene in as.vector(annot$LOCUS)) {
# 		gene_len_dict[[gene]] <- annot[annot$LOCUS==gene,]$LENGTH
# 	}
# 	## fill gene length list in corresponding to the count matrix
# 	genes_len <- c()
# 	for (i in seq(length(genes))) {
# 		genes_len <- c(genes_len, gene_len_dict[[genes[i]]])
# 	}
# 	return(genes_len)
# }


## main
parsed <- parse_arguments()
cat('... Loading count matrix\n')
count_mtx <- read_count_matrix(parsed$file_count_mtx)
samples <- colnames(count_mtx)
genes <- rownames(count_mtx)
cat('... Noramlizing count matrix\n')
norm_count_mtx <- normalize_count(count_mtx, samples)
cat('... Writing output\n')
write_count_matrix_file(parsed$output, norm_count_mtx)
# write_indiv_count_file(parsed$output, norm_count_mtx, samples)
