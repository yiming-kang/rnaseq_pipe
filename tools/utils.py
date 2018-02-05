#!/usr/bin/python
import os
import re
import numpy as np
from itertools import combinations
import yaml


def decompose_status2bit(n):
	"""
	Decompose the bit status 
	"""
	if n == 0:
		return None
	decomp = []
	for i in np.arange(np.floor(np.log2(n)),-1,-1):
		if n - 2**i >= 0:
			decomp.append(i)
			n -= 2**i
	return decomp


def make_combinations(lst):
	"""
	Make all possible replicate combinations
	"""
	if len(lst) < 2:
		return [lst]
	combo = []
	for i in range(len(lst), 1, -1):
		for s in combinations(lst, i):
			combo.append(s)
	return combo


def load_config(json_file):
	"""
	Load configuration file (JSON) for QC thresholding and scoring
	"""
	with open(json_file) as json_data:
		d = yaml.safe_load(json_data)
	return d


def parse_gtf(filename):
	"""
	Convert gtf into dictionary for filling bed fields
	"""
	bed_dict = {}
	reader = open(filename, 'r')
	for line in reader.readlines():
		if not line.startswith('#'):
			## get info from each line
			line_split = line.split("\t")
			chrm, ltype, strand, annot = line_split[0], line_split[2], line_split[6], line_split[8]
			coords = [int(x) for x in line_split[3:5]]
			gene_id = re.search(r'gene_id "(.*)"; transcript_id "(.*)"', annot).group(1)
			## fill dictionary
			if gene_id not in bed_dict.keys():
				bed_dict[gene_id] = {'chrm': chrm, 'strand': strand, 'coords':[0,0]}
			if ltype in 'start_codon':
				if strand == '+':
					bed_dict[gene_id]['coords'][0] = min(coords)
				elif strand == '-':
					bed_dict[gene_id]['coords'][1] = max(coords)
			elif ltype in 'stop_codon':
				if strand == '+':
					bed_dict[gene_id]['coords'][1] = max(coords)
				elif strand == '-':
					bed_dict[gene_id]['coords'][0] = min(coords)
	reader.close()
	return bed_dict


def parse_gff3():
	"""
	Convert gff3 into dictionary for filling bed fields
	"""
	bed_dict = {}
	reader = open(filename, 'r')
	for line in reader.readlines():
		if not line.startswith('#'):
			## get info from each line
			line_split = line.split("\t")
			chrm, ltype, strand, annot = line_split[0], line_split[2], line_split[6], line_split[8]
			coords = [int(x) for x in line_split[3:5]]
			if ltype in 'gene':
				gene_id = re.search(r'ID=(.*);Name=(.*)', annot).group(1)
				## fill dictionary
				bed_dict[gene_id] = {'chrm': chrm, 'strand': strand, 'coords':[min(coords), max(coords)]}
	reader.close()
	return bed_dict


def mkdir_p(d):
	if not os.path.exists(d):
		os.makedirs(d)


def check_dir(d):
	if not d.endswith('/'):
		d += '/'
	return d
	