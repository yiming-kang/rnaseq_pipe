# crypto_pipe.dev

![pipeline](pipeline_illustration.png)

### REQUIREMENT

This pipeline uses SLURM workload manager to streamline the RNAseq analysis. The following tools/modules are required. Tested on the respective versions. 
	
	* SLURM v17.02.6
	* novoalign v3.07.00
	* stringtie v1.3.3b  
	* samtools v1.6
	* pandas v0.20.3
	* pysam v0.11.0
	* igv v2.3.60
	* java

### SETUP
	
1. Build index for the reference genome. Make sure the version of aligner used for genome indexing is consistent with that for read alignment. Check tool manual for details.
	
	```
	ml novoalign/3.07.00
	novoindex H99/crNeoH99.nix H99/crNeoH99.fasta 
	```

2. Make direcotries. 

	```
	mkdir -p {alignment/{novoalign},expression/{stringtie,stringtie_fpkm},job_scripts,log,reports,sequence}
	```

3. [Optional] Generate and configure the IGV genome file of the species of interest for automated IGV snapshot. 

	1. On your local computer, make a directory `$HOME/igv/<strain>`, and put in genome sequence (`.fasta`) and gene annotation (`.gtf/gff`).
	2. Open IGV app, go to Genomes > Create .genome File, load the files as instructed, and save output at `$HOME/igv/genomes/`.
	3. Copy your locally created genome file and `user-defined-genomes.txt` file at `$HOME/igv/genomes/` to the server directory `$HOME/igv/genomes/`. 
	4. Copy your local directory `$HOME/igv/<strain>` (including an indexing file `.fasta.fai` built by IGV) to the cluster directory `$HOME/igv/<strain>`.
	5. Go to the cluster IGV directory. Edit the following lines in file `$HOME/igv/prefs.properties`:

	```
	DEFINE_GENOME_INPUT_DIRECTORY_KEY=<your_cluster_home_dir>/igv/<strain>
	DEFAULT_GENOME_KEY=<strain>
	```

4. [Optional] Install NOISeq package from Bioconductor.

	```
	source("https://bioconductor.org/biocLite.R")
	biocLite("NOISeq")
	```

### USAGE

1. Preparation of sequence files and metadata file 
	
	1. Make soft link or copy fastq files to sequence. `*.gz` sequence files are acceptable.
	2. Update sample summary metadata file.

		Single batch:
		```
		python tools/prepare_samples.py -m metadata/EXPERIMENT_9.xlsx
		```

		Mutliple batches to be analyzed together:
		```
		python tools/prepare_samples.py -m metadata/EXPERIMENT_9.xlsx,metadata/EXPERIMENT_10.xlsx,metadata/EXPERIMENT_11.xlsx -g 1
		```

2. Reads alignment and expression quantification
	
	This builds the SLURM job script from sample summary. Each job requires 8 CPUs and 24GB of memory. It allows 32 jobs at maximum running in parallel, depending on the available resources (e.g. CPUs and memories). The system may also send notification to user when the run fails or completes.
	
	```
	python tools/build_stage1.py -i H99/crNeoH99.nix -r H99/crNeoH99.gtf -l H99/gids -g 10 -o job_scripts/stage1.sbatch
	sbatch job_scripts/stage1.sbatch
	```

3. Quality assessment

	1. Assess the quality of each sample. The metrics used in this assessment are:
		* Total read count
		* Percentage of uniquely aligned reads
		* Efficiency of gene perturbation
		* Replicate concordance
		* Efficiency of the replaced drug-marker gene. 

		The status is in bit form (encoding of the corresponding flags is stored in qc_config.yaml).
	
	```
	ml pandas/0.20.3
	python tools/assess_quality.py -l H99/gids -g 10 -r 4 -w CNAG_00000 -c CNAG_G418,CNAG_NAT -o reports/sample_quality.group_10.xlsx
	```

	2. [Optional] Assess the efficiency of gene perturbation. Make automated IGV snapshot of the problematic mutant and marker genes. The output snapshot is titled as `[sample]gene_mutant.png`.

	```
	ml pysam/0.11.0
	python tools/build_igv_snapshot.py -q reports/sample_quality.group_10.xlsx -g H99/crNeoH99.gtf -gm H99 -o reports/inefficient_mutation.group_10/
	sbatch job_scripts/igv_snapshot.sbatch
	```

	3. [Optional] Make saturation plots for samples grouped by genotype. The output plot is title as `[genotype].png`

	```
	ml R/3.2.1
	Rscript tools/make_saturation_curve.r -i expression/stringtie_gene_count_matrix.csv -o reports/saturation_curves.group_9_10
	```

4. Differential expression  

