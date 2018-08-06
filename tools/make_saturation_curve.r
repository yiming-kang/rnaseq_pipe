suppressMessages(library(optparse))
suppressMessages(library(NOISeq))
options(bitmapType='cairo') ## enable X11


parse_arguments <- function(){
	option_list <- list(
		make_option(c('-i', '--file_count_mtx'), 
					help='Read count matrix (genes x samples)'),
		make_option(c('-k', '--detection_k',
					help='Features detected with > k counts')),
		make_option(c('-o', '--dir_output'), 
					help='Table of DE genes ranked by ajusted p-value'))
	opt <- parse_args(OptionParser(option_list=option_list))
	return(opt)
}


load_count_matrix <- function(filepath) {
	## Load count matrix
	cnts <- read.table(filepath, sep=',', row.names=1, header=TRUE, check.names=FALSE)
	return(cnts)
}


get_genotypes <- function(cnts) {
	## Get a hash for samples to each genotype
	samples <- colnames(cnts)
	geno_hash <- list()
	for (sample in samples) {
		geno <- unlist(strsplit(sample, '-'))[1]
		if (!(geno %in% names(geno_hash))) geno_hash[[geno]] <- c()
		geno_hash[[geno]] <- c(geno_hash[[geno]], sample)
	} 
	return(geno_hash)
}


make_saturation_curves <- function(cnts, geno_hash, dir_output, det_k=0) {
	## Make satruation curves plot for each genotype
	print(names(geno_hash)) 
	for (geno in names(geno_hash)) {
		## iterate thru genotypes
		cat('... working on ', geno, '\n')
		samples <- geno_hash[[geno]]
		num_samples <- length(samples)
		## process count data
		factors <- data.frame(reps=seq(num_samples), row.names=samples)
		data <- readData(data=cnts[samples], factors=factors)
		## estimate saturation
		sat_est <- dat(data, k=det_k, ndepth=10, type='saturation')
		## plot curves
		filepath <- paste0(dir_output, '/', geno, '.png')
		png(filepath, width=960, height=960, res=150)
		explo.plot(sat_est, samples=1:num_samples)
		dev.off()
	}
}

## main
parsed_opt <- parse_arguments()
cnts <- load_count_matrix(parsed_opt$file_count_mtx)
geno_hash <- get_genotypes(cnts)
dir.create(parsed_opt$dir_output, showWarnings=FALSE)
make_saturation_curves(cnts, geno_hash, parsed_opt$dir_output, parsed_opt$detection_k)


