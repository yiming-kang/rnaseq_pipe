# crypto_pipe.dev

### REQUIREMENT
1. This pipeline uses SLURM workload manager to streamline the RNAseq analysis. Tested on SLURM v17.02.6.

2. The following tools/modules are required. Tested on the respective versions. Install if running on other clusters. 

	```
	module load novoalign/3.07.00
	module load stringtie/1.3.3b  
	module load samtools/1.6
	module load pandas/0.20.3
	module load pysam/0.11.0
	module load igv/2.3.60
	module load java
	```

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

3. [Optional] IGV-snapshot-automator needs to be installed to make IGV batch file for automated snapshot. 

	1. Clone the following repo and run the testing demo to make sure the batch file (.bat) and .png images can be generated.

	```
	cd tools/
	git clone https://github.com/stevekm/IGV-snapshot-automator.git
	cd IGV-snapshot-automator/
	bash ../interactive.sh
	module load java
	python make_IGV_snapshots.py test_data/test_alignments.bam test_data/test_alignments2.bam -bin /opt/apps/igv/2.3.60/igv.jar
	```

	2. Generate and configure the IGV genome file of the species of interest. 

		1. On your local computer, make a directory `$HOME/igv/<strain>`, and put in genome sequence (`.fasta`) and gene annotation (`.gtf/gff`).
		2. Open IGV app, go to Genomes > Create .genome File, load the files as instructed, and save output at `$HOME/igv/genomes/`.
		3. Copy your locally created genome file and `user-defined-genomes.txt` file at `$HOME/igv/genomes/` to the server directory `$HOME/igv/genomes/`. 
		4. Copy your local directory `$HOME/igv/<strain>` (including an indexing file `.fasta.fai` built by IGV) to the cluster directory `$HOME/igv/<strain>`.
		5. Go to the cluster IGV directory. Edit the following lines in file `$HOME/igv/prefs.properties`:

		```
		DEFINE_GENOME_INPUT_DIRECTORY_KEY=<your_cluster_home_dir>/igv/<strain>
		DEFAULT_GENOME_KEY=<strain>
		```


### USAGE

1. Preparation of sequence files and metadata file 
	
	1. Make soft link or copy fastq files to sequence.
	2. Update sample summary metadata file.

2. Reads alignment and expression quantification
	
	This builds the SLURM job script from sample metadata. Each job requires 8 CPUs and 24GB of memory. It allows 32 jobs at maximum running in parallel, depending on the available resources (e.g. CPUs and memories). The system may also send notification to user when the run fails or completes.
	
	```
	python tools/build_stage1.py -s metadata/sample_summary.txt -i H99/crNeoH99.nix -r H99/crNeoH99.gtf -l H99/gids -g 10 > job_scripts/stage1.sbatch
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
	ml pysam/0.11.0
	python tools/assess_quality.py -s metadata/sample_summary.txt -l H99/gids -g 10 -w CNAG_00000 -c CNAG_G418,CNAG_NAT -o reports/sample_quality.group_10.txt
	```

	2. [Optional] Assess the efficiency of gene perturbation. Make automated IGV snapshot of the problematic mutant and marker genes. The output snapshot is titled as `[sample]gene_mutant.png`.

	```
	python tools/build_igv_snapshot.py -q reports/sample_quality.group_10.txt -g H99/crNeoH99.gtf -gm H99 -o reports/inefficient_mutation.group_10/
	sbatch job_scripts/igv_snapshot.sbatch
	```

4. Differential expression  

