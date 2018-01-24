# crypto_pipe.dev

### Usage
1. Stage 1: reads alignment, gene expression quantification

	1. Configure job script 
	```
	python tools/build_stage1.py -s metadata/sample_summary.txt \
	-i H99/crNeo99.nix -r H99/crNeo99.gtf -g H99/gids \
	> job_scripts/stage_1.makefile 
	```
	2. Submit SLURM job
	```
	sbatch job_scripts/run_stage_1.sbatch
	```

