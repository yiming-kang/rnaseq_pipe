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
	parser.add_argument('-o', '--output_dir', required=True,
						help='Directory for inefficient mutant samples.')
	parser.add_argument('--qc_configure', default='tools/qc_config.yaml',
						help='Configuration file for quality assessment.')
	return parser.parse_args(argv[1:])


def find_inefficient_mutation(df):
	"""
	Find the samples and gene that are incompletely perturbed
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
				ineffmut_dict[sample] = []
				for j in range(len(mutants)):
					if mutants[j].endswith('over'): 
						if mutfows[j] < QC_dict['over_fow']['threshold']:
							ineffmut_dict[sample].append(mutants[j])
					else:
						if mutfows[j] > QC_dict['del_fow']['threshold']:
							ineffmut_dict[sample].append(mutants[j])
	return ineffmut_dict


def index_bams(in_dict, aligner='novoalign'):
	"""
	Index the bam files in the same directory of bam
	"""
	for sample in in_dict.keys():
		bam_file = '/'.join(['alignment', aligner, sample, 
							'aligned_reads_sorted.bam']) 
		print '... indexing', bam_file
		out = pysam.index(bam_file)


def create_fields(in_dict, flank=500):
	"""
	"""
	pass


def write_file(in_dict, out_fn):
	writer = open(out_fn, 'w')
	for k in sorted(in_dict):
		writer.write('%s\t%s\n' % (k, '.'.join(in_dict[k])))
	writer.close()


def main(argv):
	parsed = parse_args(argv)

	global QC_dict
	QC_dict = load_config(parsed.qc_configure)
	df = pd.read_csv(parsed.sample_quality, delimiter='\t', dtype=np.str)
	
	if not os.path.exists(parsed.output_dir):
		os.makedirs(parsed.output_dir)
	inefficient_mutant_dict = find_inefficient_mutation(df)
	index_bams(inefficient_mutant_dict)

	output_filename = parsed.output_dir +'/samples_mutants.txt'	
	write_file(inefficient_mutant_dict, output_filename)
	


if __name__ == '__main__':
	main(sys.argv)