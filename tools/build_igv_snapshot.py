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
	parser.add_argument('-a', '--gene_annotation', required=True,
						help='Gene annotation file in GTF/GFF.')
	parser.add_argument('-g', '--group_num', required=True,
						help='Experiment group number.')
	parser.add_argument('-gm', '--igv_genome', required=True,
						help='Genome created in IGV.')
	parser.add_argument('-w', '--wildtype', required=True,
						help='Wildtype name.')
	parser.add_argument('-o', '--output_dir', required=True,
						help='Directory for inefficient mutant samples.')
	parser.add_argument('--qc_configure', default='tools/qc_config.yaml',
						help='Configuration file for quality assessment.')
	parser.add_argument('--igv_flank', default=500, type=int,
						help='Flanking region to add to the region of interest in IGV snapshot.')
	parser.add_argument('--snapshot_format', default='png',
						help='Image format of IGV snapshot. Choose from png (deafult) or svg. The latter format offers vector image.')
	parser.add_argument('--mail_user',
						help='Email address to send notifications on the jobs.')
	return parser.parse_args(argv[1:])


def find_inefficient_mutation(df, wt):
	"""
	Find the samples and genes that are incompletely perturbed
	"""
	ineffmut_dict = {}
	ineffmuts = set()
	## add inefficient mutants to mutant samples
	for i,row in df.iterrows():
		sample = row['GENOTYPE'] +'-'+ str(row['SAMPLE'])
		mutants = row['GENOTYPE'].split('.')
		if type(row['MUT_FOW']) == np.str:
			mutfows = [float(x) for x in row['MUT_FOW'].split(',')]
			## find the corresponding bit
			bit_comp = decompose_status2bit(int(row["STATUS"]))
			if bit_comp and np.log2(QC_dict['MUT_FOW']['DELETION']['status']) in bit_comp:
				## find which gene is inefficiently perturbed
				ineffmut_dict[sample] = {'gene':[], 'bam':None, 'bed':None}
				for j in range(len(mutants)):
					if mutants[j].endswith('over'): 
						if mutfows[j] < QC_dict['MUT_FOW']['OVEREXPRESSION']['threshold']:
							ineffmut_dict[sample]['gene'].append(mutants[j])
							ineffmuts.add(mutants[j])
					else:
						if mutfows[j] > QC_dict['MUT_FOW']['DELETION']['threshold']:
							ineffmut_dict[sample]['gene'].append(mutants[j])
							ineffmuts.add(mutants[j])
	## add inefficient mutants to wildtype samples
	for i,row in df.iterrows():
		if row['GENOTYPE'] == wt:
			sample = row['GENOTYPE'] +'-'+ str(row['SAMPLE'])
			ineffmut_dict[sample] = {'gene':list(ineffmuts), 'bam':None, 'bed':None}
	return ineffmut_dict


def index_bams(ineffmut_dict, aligner='novoalign'):
	"""
	Index the bam files in the same directory of bam, and add bam filepath
	"""
	for sample in ineffmut_dict.keys():
		bam_file = '/'.join(['alignment', aligner, sample, 'aligned_reads_sorted.bam'])
		ineffmut_dict[sample]['bam'] = bam_file
		if not os.path.isfile(bam_file +'.bai'):
			print '\t Indexing', bam_file
			out = pysam.index(bam_file)
	return ineffmut_dict


def create_igv_region(ineffmut_dict, gene_annot, igv_output_dir, flank, fig_format='png'):
	"""
	Create bed files to describe IGV region of interest
	"""
	## get gene dictionary with chromsom, gene coordinates, strand
	if gene_annot.endswith('gtf'):
		gene_annot_dict = parse_gtf(gene_annot)
	elif gene_annot.endswith('gff') or gene_annot.endswith('gff3'):
		gene_annot_dict = parse_gff3(gene_annot)
	else:
		sys.exit("ERROR: The gene annotation format cannot be recognized.")
	## create gene body region bed file
	for sample in ineffmut_dict.keys():
		igv_bed_filepath = igv_output_dir + sample +'.bed'
		ineffmut_dict[sample]['bed'] = igv_bed_filepath
		writer = open(igv_bed_filepath, 'w')
		for gene in ineffmut_dict[sample]['gene']:
			d = gene_annot_dict[gene]
			writer.write('%s\t%d\t%d\t[%s]%s.%s\n' % \
				(d['chrm'], d['coords'][0]-flank, d['coords'][1]+flank, sample, gene, fig_format))
	return ineffmut_dict


def write_job_script(ineffmut_dict, igv_genome, igv_output_dir, fig_format='png', email=None, job_script='job_scripts/igv_snapshot.sbatch'):
	"""
	Write job script to make IGV snapshot
	"""
	num_samples = len(ineffmut_dict.keys())
	job = '#!/bin/bash\n#SBATCH -N 1\n#SBATCH --mem=5G\n'
	job += '#SBATCH -D ./\n#SBATCH -o log/igv_snapshot_%A.out\n#SBATCH -e log/igv_snapshot_%A.err\n#SBATCH -J igv_snapshot\n'
	if email is not None:
		job += '#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=%s\n' % email
	job += '\nml java\n'
	for sample in ineffmut_dict.keys():
		bam_file = ineffmut_dict[sample]['bam']
		bed_file = ineffmut_dict[sample]['bed']
		job += 'python tools/make_IGV_snapshots.py %s -bin /opt/apps/igv/2.4.7/igv.jar -nf4 -r %s -g %s -fig_format %s -o %s\n' % (bam_file, bed_file, igv_genome, fig_format, igv_output_dir)
	writer = open(job_script, 'w')
	writer.write('%s' % job)
	writer.close()


def main(argv):
	parsed = parse_args(argv)
	output_dir = check_dir(parsed.output_dir)
	mkdir_p(output_dir)

	global QC_dict
	QC_dict = load_config(parsed.qc_configure)

	print '... Parsing QA file'
	df = pd.read_excel(parsed.sample_quality, dtype=np.str)
	ineffmut_dict = find_inefficient_mutation(df, parsed.wildtype)
	print '... Checking bam indexing'
	ineffmut_dict = index_bams(ineffmut_dict)
	print '... Creating IGV snapshot scripts'
	ineffmut_dict = create_igv_region(ineffmut_dict, parsed.gene_annotation, output_dir, flank=parsed.igv_flank, fig_format=parsed.snapshot_format)
	job_filename = 'job_scripts/igv_snapshot.group_' + parsed.group_num + '.sbatch'
	write_job_script(ineffmut_dict, parsed.igv_genome, output_dir, fig_format=parsed.snapshot_format, email=parsed.mail_user, job_script=job_filename)


if __name__ == '__main__':
	main(sys.argv)