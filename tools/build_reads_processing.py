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
	parser.add_argument('--stranded', default='no',
						help='Option for strand-specific protocol. Use yes/no/reverse.')
	parser.add_argument('--feature_types', default='CDS',
						help='List of annotated features (e.g. CDS and exon), delimited by ",". Default uses only CDS.')
	parser.add_argument('--cpus', default=8, type=int,
						help='CPUs per task. Default is 8.')
	parser.add_argument('--mem', default=24, type=int,
						help='Memory usage in G. Default is 24.')
	parser.add_argument('--mail_user',
						help='Email address to send notifications on the jobs.')
	return parser.parse_args(argv[1:])


def build_header(samples_filepath, group, primary_feature="CDS", cpus=8, mem=24, email=None):
	"""
	Build scripts for SBATCH configurations
	"""
	## parse sample sheet
	samples = pd.read_excel(samples_filepath, dtype=np.str)
	samples_valid = samples[(samples['GROUP'] == group) & (samples['ST_PIPE'] != '1')]
	n_lines = samples_valid.shape[0]
	## write a lookup file
	lookup_filepath = prepare_lookup_file(samples_valid, group, primary_feature)
	## write job script
	job = '#!/bin/bash\n#SBATCH -N 1\n#SBATCH --cpus-per-task=%d\n#SBATCH --mem=%dG\n' % (cpus, mem)
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


def build_expression_quantification(reference_gtf, feature_types, expr_tool='htseq', stranded='no'):
	"""
	Build job scripts for quantifying gene expression levels using StringTie or HTSeq
	"""
	is_gff = True if reference_gtf.split('.')[-1] == 'gff' else False
	job = 'if [[ ! -z ${data1} && ! -z ${data2} ]]; then\n' \
			+ '\tmkdir -p expression/htseq/${data1}\n'
	if expr_tool == 'htseq':
		for feature_type in feature_types:
			if is_gff:
				job += '\thtseq-count -f bam -i ID -s %s -t %s alignment/novoalign/${data1}/aligned_reads_sorted.bam %s > expression/htseq/${data1}/%s_count.tsv\n' % (stranded, feature_type, reference_gtf, feature_type.lower())
			else:
				job += '\thtseq-count -f bam -s %s -t %s alignment/novoalign/${data1}/aligned_reads_sorted.bam %s > expression/htseq/${data1}/%s_count.tsv\n' % (stranded, feature_type, reference_gtf, feature_type.lower())
	elif expr_tool == 'stringtie':		
		job += '\tstringtie -p ${SLURM_CPUS_PER_TASK} alignment/novoalign/${data1}/aligned_reads_sorted.bam -G %s -e -o expression/stringtie/${data1}/stringtie_out.gtf -A expression/stringtie/${data1}/gene_abundances.tab\n' % reference_gtf
	job += 'fi\n\n'
	return job


def prepare_lookup_file(samples, group, primary_feature, expr_tool='htseq'):
	"""
	Prepare lookup file for sbatch runs
	"""
	expr_tool_dict = {'htseq': '%s_count.tsv' % primary_feature,
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
	reference_gtf_suffix = parsed.reference_gtf.split('.')[-1]
	if reference_gtf_suffix not in ['gtf', 'gff']:
		sys.exit('ERROR: wrong gene annotation format.')
	if not os.path.exists(parsed.reference_gtf):
		sys.exit('ERROR: %s does not exist.' % parsed.reference_gtf)
	if not os.path.exists(os.path.dirname(parsed.output_filepath)):
		sys.exit('ERROR: %s does not exist.' % os.path.dirname(parsed.output_filepath))
	if not parsed.stranded in ['yes', 'no', 'reverse']:
		sys.exit('ERROR: wrong value for stranded argument.')

	feature_types = [x.strip() for x in parsed.feature_types.split(",")]
	
	print '... Building header'
	jobs = build_header(parsed.samples, parsed.group_num, feature_types[0],
						cpus=parsed.cpus, mem=parsed.mem, email=parsed.mail_user)
	print '... Building scripts for reads alignment'
	jobs += build_alignment(parsed.genome_index)
	print '... Building scripts for gene expression quantification'
	jobs += build_expression_quantification(parsed.reference_gtf, feature_types, 
						stranded=parsed.stranded)
	write_file(jobs, parsed.output_filepath)
	

if __name__ == '__main__':
	main(sys.argv)
