
suppressMessages(library(xlsx))

import_deseq2 <- function() suppressMessages(library(DESeq2))
import_edger <- function() suppressMessages(library(edgeR))


parse_metadata <- function(design_filepath, qa_filepath) {
	## Parse design table and sample quality sheet, and output contrast dictionary
	## load files
	design <- read.xlsx(design_filepath, sheetIndex=1, check.names=FALSE, na.string="NA")
	qa <- read.xlsx(qa_filepath, sheetIndex=1, check.names=FALSE)
	## find valid samples
	valid_samples <- c()
	for (i in 1:nrow(qa)) {
		if (!is.na(qa$MANUAL_AUDIT[i])) { 
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
			if (is.na(design[i,col]))
				next
			if (design[i,col] == 0 & sample_id %in% valid_samples)
				contrast_dict[[col]][['0']] <- c(contrast_dict[[col]][['0']], paste(design[i,'GENOTYPE'], sample_id, sep='-'))
			else if (design[i,col] == 1 & sample_id %in% valid_samples)
				contrast_dict[[col]][['1']] <- c(contrast_dict[[col]][['1']], paste(design[i,'GENOTYPE'], sample_id, sep='-'))
		}
		## remove disqualified contrast group, where no replicate is available for both samples
		if(length(contrast_dict[[col]][['0']]) < 2 & length(contrast_dict[[col]][['1']]) < 2) {
			contrast_dict[[col]] <- NULL
			cat('WARNING: Contrast group', col, 'will not be used, due to no replicate for both samples.\n')
		}
	}
	return(contrast_dict)
}


run_deseq2 <- function(cnt_mtx, contrast_dict, header, output_dir) {
	## Run DESeq2 analysis of each contrast group
	## prepare count and design matrices 
	contrast <- contrast_dict[[header]]
	samples_0 <- contrast[['0']]
	samples_1 <- contrast[['1']]
	if (length(samples_0) > 0 & length(samples_1) > 0) {
		condition <- c(rep('0',length(samples_0)),rep('1',length(samples_1)))
		samples <- c(contrast[['0']], contrast[['1']])
		coldata <- data.frame(condition, row.names=samples)
		## prepare DESeq2 dataset
		dds <- DESeqDataSetFromMatrix(countData=cnt_mtx[samples], 
									colData=coldata, design=~condition)
		## TODO: filter low count genes (as in EdgeR)
		## normalize with in contrast group
		dds <- estimateSizeFactors(dds)
		## DESeq2 testing
		dds <- DESeq(dds)
		res <- results(dds)
		res <- res[order(res$padj),]
		## write result
		filepath <- paste0(output_dir, '/', header, '.txt')
		write.table(res, file=filepath, quote=FALSE, sep='\t')
	}
}


run_edger <- function(cnt_mtx, contrast_dict, header, output_dir, mode='classic') {
	## Run EdgeR analysis of each contrast group
	## prepare count and design matrices
	contrast <- contrast_dict[[header]]
	samples_0 <- contrast[['0']]
	samples_1 <- contrast[['1']]
	if (length(samples_0) > 0 & length(samples_1) > 0) {
		## prepare EdgeR dataset
		condition <- c(rep('0',length(samples_0)),rep('1',length(samples_1)))
		samples <- c(contrast[['0']], contrast[['1']])
		dds <- DGEList(counts=cnt_mtx[samples], group=factor(condition))
		design_mtx <- model.matrix(~0 + dds$samples$group)
		# design_mtx <- model.matrix(~dds$samples$group)
		colnames(design_mtx) <- levels(dds$samples$group)
		## filter low count genes
		## valid genes should have CPM > 1 in at least two samples per group
		# valid_genes <- rowSums(cpm(dds) > 1) >= 2
		valid_genes <- rowSums(cpm(dds) > 2) >= length(samples)
		dds <- dds[valid_genes, keep.lib.sizes=FALSE]
		## estimate lib size and norma factor
		dds <- calcNormFactors(dds)
		if (mode == 'classic') {
			## estimate dispersion
			dds <- estimateDisp(dds, design_mtx)
			## exact test
			lrt <- exactTest(dds)
		} else if (mode == 'complex') {
			## estimate GLM dispersion
			dds <- estimateGLMCommonDisp(dds, design_mtx)
			dds2 <- estimateGLMTrendedDisp(dds, design_mtx, method='auto')
			dds2 <- estimateGLMTagwiseDisp(dds2, design_mtx)
			## GLM testing for DE
			fit <- glmFit(dds2, design_mtx)
			lrt <- glmLRT(fit)
			# fit <- glmQLFit(dds, design_mtx)
			# lrt <- glmQLFTest(fit)
		}
		res <- topTags(lrt, n=dim(cnt_mtx)[1], adjust.method='BH', sort.by='PValue')
		## write result
		filepath <- paste0(output_dir, '/', header, '.txt')
		write.table(res, file=filepath, quote=FALSE, sep='\t')
	}	
}
