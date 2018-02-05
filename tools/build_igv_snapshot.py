#!/usr/bin/python
import sys
import os
import argparse
import pandas as pd
import numpy as np
import pysam
from utils import *


def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-q', '--sample_quality', required=True,
						help='Sample quality file.')
	parser.add_argument('-g', '--gene_annotation', required=True,
						help='Gene annotation file in GTF/GFF.')
	parser.add_argument('-gm', '--igv_genome', required=True,
						help='Genome created in IGV.')
	parser.add_argument('-o', '--output_dir', required=True,
						help='Directory for inefficient mutant samples.')
	parser.add_argument('--qc_configure', default='tools/qc_config.yaml',
						help='Configuration file for quality assessment.')
	parser.add_argument('--igv_flank', default=500, type=int,
						help='Flanking region to add to the region of interest in IGV snapshot.')
	parser.add_argument('--mail_user',
						help='Email address to send notifications on the jobs.')
	return parser.parse_args(argv[1:])


def find_inefficient_mutation(df):
	"""
	Find the samples and genes that are incompletely perturbed
	"""
	ineffmut_dict = {}
	for i,row in df.iterrows():
		sample = row['GENOTYPE'] +'-'+ str(row['SAMPLE'])
		mutants = row['GENOTYPE'].split('.')
		if type(row['MUT_FOW']) == np.str:
			mutfows = [float(x) for x in row['MUT_FOW'].split(',')]
			## find the corresponding bit
			bit_comp = decompose_status2bit(int(row["STATUS"]))
			if bit_comp and np.log2(QC_dict['del_fow']['status']) in bit_comp:
				## find which gene is inefficiently perturbed
				ineffmut_dict[sample] = {'gene':[], 'bam':None, 'bed':None}
				for j in range(len(mutants)):
					if mutants[j].endswith('over'): 
						if mutfows[j] < QC_dict['over_fow']['threshold']:
							ineffmut_dict[sample]['gene'].append(mutants[j])
					else:
						if mutfows[j] > QC_dict['del_fow']['threshold']:
							ineffmut_dict[sample]['gene'].append(mutants[j])
	return ineffmut_dict


def index_bams(ineffmut_dict, aligner='novoalign'):
	"""
	Index the bam files in the same directory of bam, and add bam filepath
	"""
	for sample in ineffmut_dict.keys():
		bam_file = '/'.join(['alignment', aligner, sample, 'aligned_reads_sorted.bam'])
		ineffmut_dict[sample]['bam'] = bam_file
		if not os.path.isfile(bam_file +'.bai'):
			print '... indexing', bam_file
			out = pysam.index(bam_file)
	return ineffmut_dict


def create_igv_region(ineffmut_dict, gene_annot, igv_output_dir, flank):
	"""
	Create bed files to describe IGV region of interest
	"""
	## get gene dictionary with chromsom, gene coordinates, strand
	if gene_annot.endswith('gtf'):
		gene_annot_dict = parse_gtf(gene_annot)
	elif gene_annot.endswith('gff'):
		gene_annot_dict = parse_gff(gene_annot)
	else:
		sys.exit("Error: Incorrect gene annotation file.")
	## create gene body region bed file
	for sample in ineffmut_dict.keys():
		igv_bed_filepath = igv_output_dir + sample +'.bed'
		ineffmut_dict[sample]['bed'] = igv_bed_filepath
		writer = open(igv_bed_filepath, 'w')
		for gene in ineffmut_dict[sample]['gene']:
			d = gene_annot_dict[gene]
			writer.write('%s\t%d\t%d\t[%s]%s.png\n' % \
				(d['chrm'], d['coords'][0]-flank, d['coords'][1]+flank, sample, gene))
	return ineffmut_dict


def write_job_script(ineffmut_dict, lookup_filepath, igv_genome, igv_output_dir, email=None, job_script='job_scripts/igv_snapshot.sbatch'):
	"""
	Write job script to make IGV snapshot
	"""
	num_samples = len(ineffmut_dict.keys())
	job = '#!/bin/bash\n#SBATCH -N 1\n#SBATCH --mem=5G\n'
	job	+= '#SBATCH --array=1-%d\n' % num_samples 
	job += '#SBATCH -D ./\n#SBATCH -o log/igv_snapshot_%A_%a.out\n#SBATCH -e log/igv_snapshot_%A_%a.err\n#SBATCH -J igv_snapshot\n'
	if email is not None:
		job += '#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=%s\n' % email
	job += '\nml java\n'
	job += 'read sample_id genes bam_file bed_file < <( sed -n ${SLURM_ARRAY_TASK_ID}p %s )\nset -e\n\n' % lookup_filepath
	job += 'python tools/IGV-snapshot-automator/make_IGV_snapshots.py $bam_file -bin /opt/apps/igv/2.4.7/igv.jar -nf4 -g %s -r $bed_file -o %s' % (igv_genome, igv_output_dir)
	writer = open(job_script, 'w')
	writer.write('%s' % job)
	writer.close()


def write_file(ineffmut_dict, out_fn):
	"""
	Write the sample and genes that are inefficiently perturbed
	"""
	writer = open(out_fn, 'w')
	for sample in sorted(ineffmut_dict):
		writer.write('%s\t%s\t%s\t%s\n' % (sample,
					'.'.join(ineffmut_dict[sample]['gene']),
					ineffmut_dict[sample]['bam'],
					ineffmut_dict[sample]['bed']))
	writer.close()


def main(argv):
	parsed = parse_args(argv)
	output_dir = check_dir(parsed.output_dir)
	mkdir_p(output_dir)

	global QC_dict
	QC_dict = load_config(parsed.qc_configure)

	df = pd.read_csv(parsed.sample_quality, delimiter='\t', dtype=np.str)
	ineffmut_dict = find_inefficient_mutation(df)
	ineffmut_dict = index_bams(ineffmut_dict)
	ineffmut_dict = create_igv_region(ineffmut_dict, parsed.gene_annotation, output_dir, parsed.igv_flank)

	output_filename = output_dir +'samples_mutants.txt'
	write_job_script(ineffmut_dict, output_filename, parsed.igv_genome, output_dir, parsed.mail_user)
	write_file(ineffmut_dict, output_filename)	


if __name__ == '__main__':
	main(sys.argv)