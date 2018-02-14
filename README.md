# crypto_pipe.dev

This RNA-seq analysis pipeline is designed for processing data generated from gene perturbation, time series, and multiple level treatment experiments. As illustrated in the following figure, the pipeline modules align Illumina reads, quantify transcriptomic expression, assess sample quality, analyze differential expression (DE), and generate user-friendly reports. Its optional modules provide more thorough quality analysis and help guide future experimental design. The pipeline aim at minimizing human effort in metadata maintenance, sample quality assessment and DE design, while maximizing the versatility in handling complex experiment design. 

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
novoindex <genome>.nix <genome>.fasta 
```

2. Make direcotries. 

```
mkdir -p {alignment/{novoalign},expression/{stringtie,stringtie_count_matrix},job_scripts,log,reports,sequence}
```

3. **[Optional]** Generate and configure the IGV genome file of the species of interest for automated IGV snapshot. 

	1. On your local computer, make a directory `$HOME/igv/<strain>`, and put in genome sequence (`.fasta`) and gene annotation (`.gtf/gff`).
	2. Open IGV app, go to Genomes > Create .genome File, load the files as instructed, and save output at `$HOME/igv/genomes/`.
	3. Copy your locally created genome file and `user-defined-genomes.txt` file at `$HOME/igv/genomes/` to the server directory `$HOME/igv/genomes/`. 
	4. Copy your local directory `$HOME/igv/<strain>` (including an indexing file `.fasta.fai` built by IGV) to the cluster directory `$HOME/igv/<strain>`.
	5. Go to the cluster IGV directory. Edit the following lines in file `$HOME/igv/prefs.properties`:

```
DEFINE_GENOME_INPUT_DIRECTORY_KEY=<your_cluster_home_dir>/igv/<strain>
DEFAULT_GENOME_KEY=<strain>
```

4. **[Optional]** Install NOISeq package from Bioconductor.

```
source("https://bioconductor.org/biocLite.R")
biocLite("NOISeq")
```

### USAGE

1. Data preparation and pre-alignment QC
	
	1. Make a subdirectory `sequence/run_<#>_samples/` for the batch. Make soft link or copy the reads files to it. Both zipped `*.gz` or unzipped fastq reads files are acceptable.
	2. **[Optional]** Update the thresholds for sample QC in configrue file `tools/qc_config.yaml`.
	3. Prepare samples that will be analyzed together in the same analysis group. This module sifts low-quality sample before alignment, update sample summary sheet, and generates a design table with default contrast group that will be used in DE analysis.

```
python tools/prepare_samples.py -g <group_#> -m metadata/<metadata>.xlsx \
		-d reports/design_table.<group_#>.xlsx -w <wildtype> 
```

2. Reads alignment and transcriptomic expression quantification
	
	This module builds the SLURM job script from sample summary. Each job requires 8 CPUs and 24GB of memory. It allows 32 jobs at maximum running in parallel, depending on the available resources (e.g. CPUs and memories). The user may opt to recive email notification when the run fails or completes.
	
```
python tools/build_stage1.py -g <group_#> -l <gene_list> -i <genome>.nix \
		-r <gene_annotation>.gtf -o job_scripts/<stage1_job>.sbatch \
		--mail_user <email_address>
sbatch job_scripts/<stage1_job>.sbatch
```

3. Quality assessment

	1. Assess the quality of each sample. The QC metrics are the followings:
		* Total read count
		* Percentage of uniquely aligned reads
		* Efficiency of gene perturbation
		* Replicate concordance
		* Efficiency of the replaced drug-marker gene. 

	The status (in bit form) summarizes the overall quality of each sample. The encoding of the corresponding metric is stored in `tools/qc_config.yaml`. `<markger_genes>` should be tab delimited, if more than one marker is used.
	
```
ml pandas/0.20.3
python tools/assess_quality.py -g <group_#> -l <gene_list> -r <max_replicate_#> \
		-w <wildtype> -c <marker_genes> -o reports/sample_quality.<group_#>.xlsx
```

	2. **[Optional]** Assess the efficiency of gene perturbation. Make automated IGV snapshot of the problematic mutant and marker genes. The output snapshot is titled `[<sample>]<gene_mutant>.png`. Two snapshots will be made for double mutant, and so forth.

```
ml pysam/0.11.0
python tools/build_igv_snapshot.py -q reports/sample_quality.<group_#>.xlsx \
		-g <gene_annotation>.gtf -gm <genome> \
		-o reports/inefficient_mutation.<group_#>/
sbatch job_scripts/igv_snapshot.sbatch
```

	3. **[Optional]** Make saturation plots for samples grouped by genotype. Each output plot titled `[genotype].png` contains all replicates/samples belonging to the same genotype.

```
ml R/3.2.1
Rscript tools/make_saturation_curve.r -i <count_matrix>.csv \
		-o reports/saturation_curves.<group_#>/
```

4. Differential expression  

