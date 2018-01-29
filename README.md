# crypto_pipe.dev

### Usage
1. Preperation.
	
	1. [One time only] Build index for the reference genome. Make sure the version of aligner used for genome indexing is consistent with that for read alignment. Check tool manual for details.
	
	```
	ml novoalign/3.07.00
	novoindex H99/crNeoH99.nix H99/crNeoH99.fasta 
	```

	2. [One time only] Make direcotries. 

	```
	mkdir -p {alignment/{novoalign},expression/{stringtie,stringtie_fpkm},job_scripts,log,reports,sequence}
	```

	3. Soft link or copy fastq files to sequence.

	4. Update sample summary metadata file.

2. Align reads to reference genome and quantify gene expression levels.
	
	This builds the SLURM job script from sample metadata. Each job requires 8 CPUs and 24GB of memory. It allows 32 jobs at maximum running in parallel, depending on the available resources (e.g. CPUs and memories). 
	
	```
	python tools/build_stage1.py -s metadata/sample_summary.txt \
		-i H99/crNeoH99.nix -r H99/crNeoH99.gtf -l H99/gids -g 10 \
		> job_scripts/stage1.sbatch
	sbatch job_scripts/stage.sbatch
	```

3. Quality assessment
	
	```
	ml pandas/0.20.3
	python tools/assess_quality.py -s metadata/sample_summary.txt -l H99/gids -g 10 -w CNAG_00000
	```

4. Analyze differential expression  

