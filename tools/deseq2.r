library(optparse)
library(xlsx)
library(DESeq2)

parse_arguments <- function(){
	option_list <- list(
		make_option(c('-c', '--count_matrix'), 
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
	## Parse design table and sample quality sheet, and output contrast dictionary
	## load files
	design <- read.xlsx(design_filepath, sheetIndex=1, check.names=FALSE, na.string="NA")
	qa <- read.xlsx(qa_filepath, sheetIndex=1, check.names=FALSE)
	## find valid samples
	valid_samples <- c()
	for (i in 1:nrow(qa)) {
		if (is.na(qa$MANUAL_AUDIT[i])) { 
			if (qa$AUTO_AUDIT[i] == 0) 
				valid_samples <- c(valid_samples, as.character(qa$SAMPLE[i]))
		} else {
			if (qa$MANUAL_AUDIT[i] == 0) 
				valid_samples <- c(valid_samples, as.character(qa$SAMPLE[i]))
		}
	}
	## group samples using constrast descriptor as key and valid samples as value
	d_cols <- colnames(design)
	contrast_dict <- list()
	for (col in d_cols) {
		## for each contrast group
		col_split <- strsplit(col, '\\[|\\]')[[1]]
		if (length(col_split) <= 1) next 
		## get the contrast type, e.g. GENOTYPE
		contrast_type <- strsplit(col_split[3], '\\:')[[1]][1]
		contrast_dict[[col]] <- list('0'=c(), '1'=c())
		for (i in 1:nrow(design)) {
			## for each sample if is assigned with a number and is valid 
			sample_id <- as.character(design[i,'SAMPLE'])
			if (design[i,col] == 0 & sample_id %in% valid_samples) 
				contrast_dict[[col]][['0']] <- c(contrast_dict[[col]][['0']], paste(design[i,contrast_type], sample_id, sep='-'))
			else if (design[i,col] == 1 & sample_id %in% valid_samples) 
				contrast_dict[[col]][['1']] <- c(contrast_dict[[col]][['1']], paste(design[i,contrast_type], sample_id, sep='-'))
		}
	}
	return(contrast_dict)
}


analyze_deseq2 <- function(cnt_mtx, contrast_dict, header, output_dir) {
	## Run DESeq2 analysis of each contrast group
	## prepare count and design matrices for DESeq2
	contrast <- contrast_dict[[header]]
	samples_0 <- contrast[['0']]
	samples_1 <- contrast[['1']]
	if (length(samples_0) > 0 & length(samples_1) > 0) {
		condition <- c(rep('0',length(samples_0)),rep('1',length(samples_1)))
		samples <- c(contrast[['0']], contrast[['1']])
		coldata <- data.frame(condition, row.names=samples)
		dds <- DESeqDataSetFromMatrix(countData=cnt_mtx[samples], 
									colData=coldata, design=~condition)
		## run DESeq2
		dds <- DESeq(dds)
		res <- results(dds)
		res <- res[order(res$padj),]
		## write result
		filepath <- paste0(output_dir, '/', header, '.txt')
		write.table(res, file=filepath, quote=FALSE, sep='\t')
	}
}

## main
parsed_opt <- parse_arguments()
dir.create(parsed_opt$output_dir, showWarnings=FALSE)
cat('... Parsing metadata\n')
contrast_dict <- parse_metadata(parsed_opt$design_table, parsed_opt$qa_table)
cat('... Loading count matrix\n')
count_matrix <- read.csv(parsed_opt$count_matrix, check.names=FALSE, row.names=1)
for (header in names(contrast_dict)) {
	cat('... Working on', header, '\n')
	analyze_deseq2(count_matrix, contrast_dict, header, parsed_opt$output_dir)
}

	