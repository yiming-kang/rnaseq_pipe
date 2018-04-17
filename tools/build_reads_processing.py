#!/usr/bin/python
import sys
import argparse
import os.path
import pandas as pd
import numpy as np


def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--samples', required=True,
						help='Sample summary metadata file.')
	parser.add_argument('-i', '--genome_index', required=True,
						help='Genome index file built for the aligner, e.g. *.nix for novoalign.')
	parser.add_argument('-r', '--reference_gtf', required=True,
						help='Annotation reference file in GTF/GFF3 format.')
	parser.add_argument('-g', '--group_num', required=True,
						help='Experiment group number.')
	parser.add_argument('-o', '--output_filepath', required=True,
						help='Filepath of sbatch job script.')
	parser.add_argument('--stranded', action='store_true',
						help='Use if it is strand-specific protocol.')
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
	n_lines = samples_valid.shape[0]
	## write a lookup file
	lookup_filepath = prepare_lookup_file(samples_valid, group)
	## write job script
	job = '#!/bin/bash\n#SBATCH -N 1\n#SBATCH --cpus-per-task=8\n#SBATCH --mem=24G\n'
	job	+= '#SBATCH --array=1-%d%%%d\n' % (n_lines, min(n_lines,32))
	job += '#SBATCH -D ./\n#SBATCH -o log/readsproc_%A_%a.out\n#SBATCH -e log/readsproc_%A_%a.err\n#SBATCH -J readsproc\n'
	if email is not None:
		job += '#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=%s\n' % email
	job += '\nml novoalign/3.07.00\nml samtools/1.6\nml stringtie/1.3.3b\nml htseq/0.9.1\nml R/3.2.1\n'
	job += 'read data1 data2 < <( sed -n ${SLURM_ARRAY_TASK_ID}p %s )\nset -e\n\n' % lookup_filepath
	return job


def build_alignment(genome_index):
	"""
	Build job scripts for aligning reads to ref genome using Novoalign
	"""
	job = 'if [[ ! -z ${data1} && ! -z ${data2} ]]; then\n'
	job += 'if [ ! -f alignment/novoalign/${data1}/aligned_reads_sorted.bam ]; then\n' \
		+ '\trm -rf alignment/novoalign/${data1} || :\n' \
		+ '\tmkdir -p alignment/novoalign/${data1}\n' 
	job += '\tnovoalign -c ${SLURM_CPUS_PER_TASK} -o SAM -d %s -f ${data2} 2>alignment/novoalign/${data1}/novoalign.log | samtools view -bS > alignment/novoalign/${data1}/aligned_reads.bam\n' % genome_index
	job += '\tnovosort --threads ${SLURM_CPUS_PER_TASK} alignment/novoalign/${data1}/aligned_reads.bam >  alignment/novoalign/${data1}/aligned_reads_sorted.bam 2> alignment/novoalign/${data1}/novosort.log\n' \
		+ '\trm alignment/novoalign/${data1}/aligned_reads.bam\n' \
		+ 'fi\n'
	job += 'fi\n\n'
	return job


def build_expression_quantification(reference_gtf, expr_tool='htseq', stranded=False):
	"""
	Build job scripts for quantifying gene expression levels using StringTie or HTSeq
	"""
	is_stranded = 'yes' if stranded else 'no'
	job = 'if [[ ! -z ${data1} && ! -z ${data2} ]]; then\n'
	if expr_tool == 'htseq':
		job += 'if [ ! -f expression/htseq/${data1}/cds_count.tsv && ! -f expression/htseq/${data1}/exon_count.tsv ]; then\n' \
			+ '\trm -rf expression/htseq/${data1} || :\n' \
			+ '\tmkdir -p expression/htseq/${data1}\n' 
		job += '\thtseq-count -f bam -s %s -t CDS alignment/novoalign/${data1}/aligned_reads_sorted.bam %s > expression/htseq/${data1}/cds_count.tsv\n' % (is_stranded, reference_gtf)
		job += '\thtseq-count -f bam -s %s -t exon alignment/novoalign/${data1}/aligned_reads_sorted.bam %s > expression/htseq/${data1}/exon_count.tsv\n' % (is_stranded, reference_gtf)
	elif expr_tool == 'stringtie':		
		job += 'if [ ! -f expression/stringtie/${data1}/stringtie_out.gtf  ]; then\n' \
			+ '\trm -rf expression/stringtie/${data1} || :\n' \
			+ '\tmkdir -p expression/stringtie/${data1}\n' 
		job += '\tstringtie -p ${SLURM_CPUS_PER_TASK} alignment/novoalign/${data1}/aligned_reads_sorted.bam -G %s -e -o expression/stringtie/${data1}/stringtie_out.gtf -A expression/stringtie/${data1}/gene_abundances.tab\n' % reference_gtf
	job += 'fi\n'
	job += 'fi\n\n'
	return job


def prepare_lookup_file(samples, group, expr_tool='htseq'):
	"""
	Prepare lookup file for sbatch runs
	"""
	expr_tool_dict = {'htseq': 'cds_count.tsv',
						'stringtie': 'stringtie_out.gtf'}
	## create lookup data for each sample
	lookup_readsproc, lookup_expr = '', ''
	for i,row in samples.iterrows(): 
		sample = row['GENOTYPE'] +'-'+ row['SAMPLE']
		expr_file = '/'.join(['expression', expr_tool, sample, expr_tool_dict[expr_tool]])
		lookup_readsproc += '%s\t%s\n' % (sample, row['FILE'])
		lookup_expr += '%s\t%s\n' % (sample, expr_file)
	## write file
	lookup_filepath_prefix = 'job_scripts/lookup_files/group_'+ group
	lookup_readsproc_filepath = lookup_filepath_prefix +'.readsproc.txt'
	lookup_expr_filepath = lookup_filepath_prefix +'.expr.txt'
	write_file(lookup_readsproc, lookup_readsproc_filepath)
	write_file(lookup_expr, lookup_expr_filepath)
	return lookup_readsproc_filepath


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
	if not os.path.exists(os.path.dirname(parsed.output_filepath)):
		sys.exit('ERROR: %s does not exist.' % os.path.dirname(parsed.output_filepath))
	
	print '... Building header'
	jobs = build_header(parsed.samples, parsed.group_num, parsed.mail_user)
	print '... Building scripts for alignment'
	jobs += build_alignment(parsed.genome_index)
	print '... Building scripts for gene expression quantification'
	jobs += build_expression_quantification(parsed.reference_gtf, stranded=parsed.stranded)
	write_file(jobs, parsed.output_filepath)
	

if __name__ == '__main__':
	main(sys.argv)
