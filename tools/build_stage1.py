#!/usr/bin/python
import sys
import argparse
import pandas as pd


def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--samples', required=True,
						help='Sample summary metadata file.')
	parser.add_argument('-i', '--genome_index', required=True,
						help='Genome index file built for the aligner, e.g. *.nix for novoalign.')
	parser.add_argument('-r', '--reference_gtf', required=True,
						help='Annotation reference file in GTF/GFF3 format.')
	parser.add_argument('-l', '--gene_list', required=True,
						help='Gene list.')
	parser.add_argument('-g', '--group_num', required=True,
						help='Experiment group number.')
	parser.add_argument('-m', '--mail_user',
						help='Email address to send notifications on the jobs.')
	return parser.parse_args(argv[1:])


def build_header(samples_filepath, group, email=None):
	## parse sample sheet
	samples = pd.read_csv(samples_filepath, delimiter='\t', dtype=str)
	samples_matched = samples[samples['GROUP'] == group]
	num_samples = samples_matched.shape[0]
	## write a lookup file
	lookup_filepath = 'job_scripts/lookup_files/group_'+ group
	writer = open(lookup_filepath, 'w')
	for i,row in samples_matched.iterrows(): 
		writer.write('%s-%s\t%s\n' % (row['GENOTYPE'], row['SAMPLE'], row['FILE']))
	writer.close()
	## write job script
	job = '#!/bin/bash\n#SBATCH -N 1\n#SBATCH --cpus-per-task=8\n#SBATCH --mem=24G\n'
	job	+= '#SBATCH --array=1-%d%%32\n' % num_samples 
	job += '#SBATCH -D ./\n#SBATCH -o log/stage1_%A_%a.out\n#SBATCH -e log/stage1_%A_%a.err\n#SBATCH -J stage1\n'
	if email is not None:
		job += '#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=%s\n' % email
	job += '\nml novoalign/3.07.00\nml samtools/1.6\nml stringtie/1.3.3b\n'
	job += 'read sample_id seq_file < <( sed -n ${SLURM_ARRAY_TASK_ID}p %s )\nset -e\n\n' % lookup_filepath
	sys.stdout.write('%s' % job)


def build_alignment(genome_index):
	job = 'if [ ! -f alignment/novoalign/${sample_id}/aligned_reads_sorted.bam ]\n' \
		+ 'then\n' \
		+ '\trm -rf alignment/novoalign/${sample_id} || :\n' \
		+ '\tmkdir -p alignment/novoalign/${sample_id}\n' 
	job += '\tnovoalign -c ${SLURM_CPUS_PER_TASK} -o SAM -d %s -f ${seq_file} 2>alignment/novoalign/${sample_id}/novoalign.log | samtools view -bS > alignment/novoalign/${sample_id}/aligned_reads.bam\n' % genome_index
	job += 'novosort --threads ${SLURM_CPUS_PER_TASK} alignment/novoalign/${sample_id}/aligned_reads.bam >  alignment/novoalign/${sample_id}/aligned_reads_sorted.bam 2> alignment/novoalign/${sample_id}/novosort.log\n' \
		+ '\trm alignment/novoalign/${sample_id}/aligned_reads.bam\n' \
		+ 'fi\n\n'
	sys.stdout.write('%s' % job)


def build_expression_quantification(reference_gtf, gene_list):
	job = 'if [ ! -f expression/stringtie_fpkm/${sample_id}.expr ]\n' \
		+ 'then\n' \
		+ '\trm -rf expression/stringtie/${sample_id} || :\n' \
		+ '\tmkdir -p expression/stringtie/${sample_id}\n' \
		+ '\trm -rf expression/stringtie_fpkm/${sample_id} || :\n' \
		+ '\tmkdir -p expression/stringtie_fpkm/${sample_id}\n' 
	job += '\tstringtie -p ${SLURM_CPUS_PER_TASK} alignment/novoalign/${sample_id}/aligned_reads_sorted.bam -G %s -e -o expression/stringtie/${sample_id}/stringtie_out.gtf -A expression/stringtie/${sample_id}/gene_abundances.tab\n' % reference_gtf
	job += '\tpython tools/stringtie2fpkm.py expression/stringtie/${sample_id}/gene_abundances.tab %s > expression/stringtie_fpkm/${sample_id}.expr\n' % gene_list
	job += 'fi\n'
	sys.stdout.write('%s' % job)


def main(argv):
	parsed = parse_args(argv)
	
	build_header(parsed.samples, parsed.group_num, parsed.mail_user)
	build_alignment(parsed.genome_index)
	build_expression_quantification(parsed.reference_gtf, parsed.gene_list)
	

if __name__ == '__main__':
	main(sys.argv)
