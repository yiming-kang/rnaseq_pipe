#!/usr/bin/python
import sys
import argparse
import os.path
import pandas as pd
import numpy as np


def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--genome_index', required=True,
						help='Genome index file built for the aligner, e.g. *.nix for novoalign.')
	parser.add_argument('-r', '--reference_gtf', required=True,
						help='Annotation reference file in GTF/GFF3 format.')
	parser.add_argument('-l', '--gene_list', required=True,
						help='Gene list.')
	parser.add_argument('-g', '--group_num', required=True,
						help='Experiment group number.')
	parser.add_argument('-o', '--output_filepath', required=True,
						help='Filepath of sbatch job script.')
	parser.add_argument('-s', '--samples', default='metadata/sample_summary.xlsx',
						help='Sample summary metadata file.')
	parser.add_argument('--read_length', default=75, type=int,
						help='Average read length.')
	parser.add_argument('--mail_user',
						help='Email address to send notifications on the jobs.')
	return parser.parse_args(argv[1:])


def build_header(samples_filepath, group, email=None):
	"""
	Build scripts for SBATCH configurations
	"""
	## parse sample sheet
	samples = pd.read_excel(samples_filepath, dtype=np.str)
	samples_valid = samples[(samples['GROUP'] == group) & (samples['ST_PIPE'] != '1')]
	n_lines = samples_valid.shape[0]+1
	## write a lookup file
	lookup_filepath = prepare_lookup_file(samples_valid, group)
	## write job script
	job = '#!/bin/bash\n#SBATCH -N 1\n#SBATCH --cpus-per-task=8\n#SBATCH --mem=24G\n'
	job	+= '#SBATCH --array=1-%d%%32\n' % n_lines
	job += '#SBATCH -D ./\n#SBATCH -o log/stage1_%A_%a.out\n#SBATCH -e log/stage1_%A_%a.err\n#SBATCH -J stage1\n'
	if email is not None:
		job += '#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=%s\n' % email
	job += '\nml novoalign/3.07.00\nml samtools/1.6\nml stringtie/1.3.3b\nml R/3.2.1\n'
	job += 'read data1 data2 < <( sed -n ${SLURM_ARRAY_TASK_ID}p %s )\nset -e\n\n' % lookup_filepath
	return job


def build_alignment(genome_index):
	"""
	Build job scripts for aligning reads to ref genome using Novoalign
	"""
	job = 'if [ ! -f alignment/novoalign/${data1}/aligned_reads_sorted.bam ]; then\n' \
		+ '\trm -rf alignment/novoalign/${data1} || :\n' \
		+ '\tmkdir -p alignment/novoalign/${data1}\n' 
	job += '\tnovoalign -c ${SLURM_CPUS_PER_TASK} -o SAM -d %s -f ${data2} 2>alignment/novoalign/${data1}/novoalign.log | samtools view -bS > alignment/novoalign/${data1}/aligned_reads.bam\n' % genome_index
	job += '\tnovosort --threads ${SLURM_CPUS_PER_TASK} alignment/novoalign/${data1}/aligned_reads.bam >  alignment/novoalign/${data1}/aligned_reads_sorted.bam 2> alignment/novoalign/${data1}/novosort.log\n' \
		+ '\trm alignment/novoalign/${data1}/aligned_reads.bam\n' \
		+ 'fi\n\n'
	return job


def build_expression_quantification(reference_gtf, gene_list):
	"""
	Build job scripts for quantifying gene expression levels using StringTie
	"""
	job = 'if [ ! -f expression/stringtie/${data1}/stringtie_out.gtf  ]; then\n' \
		+ '\trm -rf expression/stringtie/${data1} || :\n' \
		+ '\tmkdir -p expression/stringtie/${data1}\n' 
	job += '\tstringtie -p ${SLURM_CPUS_PER_TASK} alignment/novoalign/${data1}/aligned_reads_sorted.bam -G %s -e -o expression/stringtie/${data1}/stringtie_out.gtf -A expression/stringtie/${data1}/gene_abundances.tab\n' % reference_gtf
	# job += '\tpython tools/stringtie2fpkm.py expression/stringtie/${data1}/gene_abundances.tab %s > expression/stringtie_fpkm/${data1}.expr\n' % gene_list
	job += 'fi\n\n'
	return job


def build_count_normalization(group, read_length):
	"""
	Build job scripts for normalizing count matrix
	"""
	lookup_file = 'job_scripts/lookup_files/group_%s.expr.txt' % group
	job = 'if [ ! -f ${data2} ]; then\n'
	job += '\twhile true; do\n'
	job += '\t\tflag=1; while read s f; do if ! [[ -e $f ]]; then flag=0; fi; done < %s\n' % lookup_file
	job += '\t\tif [[ $flag == 1 ]]; then python tools/prepDE.py -i %s -g ${data1} -l %d; Rscript tools/normalize_counts.r -i ${data1} -o ${data2}; break; fi\n' % (lookup_file, read_length)
	job += '\t\tsleep 10\n\tdone\nfi\n'
	return job


def prepare_lookup_file(samples, group):
	"""
	Prepare lookup file for sbatch runs
	"""
	## create lookup data for each sample
	lookup_stage1, lookup_expr = '', ''
	for i,row in samples.iterrows(): 
		sample = row['GENOTYPE'] +'-'+ row['SAMPLE']
		expr_file = 'expression/stringtie/'+ sample +'/stringtie_out.gtf'
		lookup_stage1 += '%s\t%s\n' % (sample, row['FILE'])
		lookup_expr += '%s\t%s\n' % (sample, expr_file)
	## create laste line 
	lookup_stage1 += 'expression/stringtie_count_matrix/count_matrix.group_%s.csv\texpression/stringtie_count_matrix/normalized_count_matrix.group_%s.csv\n' % (group, group)
	## write file
	lookup_filepath_prefix = 'job_scripts/lookup_files/group_'+ group
	lookup_stage1_filepath = lookup_filepath_prefix +'.stage1.txt'
	lookup_expr_filepath = lookup_filepath_prefix +'.expr.txt'
	write_file(lookup_stage1, lookup_stage1_filepath)
	write_file(lookup_expr, lookup_expr_filepath)
	return lookup_stage1_filepath


def write_file(jobs, filepath):
	"""
	Write SBATCH file
	"""
	writer = open(filepath, 'w')
	writer.write('%s' % jobs)
	writer.close()


def main(argv):
	parsed = parse_args(argv)
	if not os.path.exists(parsed.samples):
		sys.exit('ERROR: %s does not exist.' % parsed.samples)
	if not os.path.exists(parsed.genome_index):
		sys.exit('ERROR: %s does not exist.' % parsed.genome_index)
	if not os.path.exists(parsed.reference_gtf):
		sys.exit('ERROR: %s does not exist.' % parsed.reference_gtf)
	if not os.path.exists(parsed.gene_list):
		sys.exit('ERROR: %s does not exist.' % parsed.gene_list)
	if not os.path.exists(os.path.dirname(parsed.output_filepath)):
		sys.exit('ERROR: %s does not exist.' % os.path.dirname(parsed.output_filepath))
	
	print '... Building header'
	jobs = build_header(parsed.samples, parsed.group_num, parsed.mail_user)
	print '... Building scripts for alignment'
	jobs += build_alignment(parsed.genome_index)
	print '... Building scripts for gene expression quantification'
	jobs += build_expression_quantification(parsed.reference_gtf, parsed.gene_list)
	print '... Building scripts for count normalization'
	jobs += build_count_normalization(parsed.group_num, parsed.read_length)
	write_file(jobs, parsed.output_filepath)
	

if __name__ == '__main__':
	main(sys.argv)
