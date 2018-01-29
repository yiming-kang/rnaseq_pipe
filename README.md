# crypto_pipe.dev

### Usage
1. Align reads to reference genome and quantify gene expression levels.

	This builds the SLURM job script from sample metadata. Each job requires 8 CPUs and 24GB of memory. It allows 32 jobs at maximum running in parallel, depending on the available resources (e.g. CPUs and memories). 

	```
	python tools/build_stage1.py -s metadata/sample_summary.txt \
		-i H99/crNeoH99.nix -r H99/crNeoH99.gtf -l H99/gids -g 10 \
		> job_scripts/stage1.sbatch
	sbatch job_scripts/stage.sbatch
	```

2. Quality assessment

	```
	ml pandas/0.20.3
	python tools/assess_quality.py -s metadata/sample_summary.txt -e 10 -l H99/gids -w CNAG_00000
	```

3. Analyze differential expression  

