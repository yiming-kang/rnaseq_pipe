#!/usr/bin/python
import sys
import argparse
import os.path
import pandas as pd
import numpy as np


def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_lookup', required=True,
						help='Lookup table of sample and the corresponding HTSeq count.')
	parser.add_argument('-o', '--output_count_file', required=True,
						help='Output count matrix.')
	parser.add_argument('-l', '--gene_list', required=True,
						help='Gene list. This gene list should be a subset or the same set of annotated genes in GTF/GFF.')
	return parser.parse_args(argv[1:])


def parse_lookup(file):
	lookup = {}
	data = np.loadtxt(file, dtype=str)
	for i in range(data.shape[0]):
		if not os.path.exists(data[i,1]):
			sys.exit('ERROR: %s does not exist.\n... Aborted preparing count matrix.' % data[i,1])
		lookup[data[i,0]] = data[i,1]
	return lookup


def create_count_matrix(lookup, gene_list):
	lookup = parse_lookup(lookup)
	gene_list = np.loadtxt(gene_list, dtype=str)
	count_mtx = gene_list.reshape(-1,1)
	header = ['gene_id']
	for sample, file in lookup.iteritems():
		print '... working on %s' % sample
		count = np.loadtxt(file, dtype=str)
		indx = np.where([gene_list[i]==count[:,0] for i in range(len(gene_list))])[1]
		count_mtx = np.hstack((count_mtx, count[indx,1].reshape(-1,1)))
		header.append(sample)
	count_mtx = np.vstack((header, count_mtx))
	return count_mtx


def main(argv):
	parsed = parse_args(argv)
	if not os.path.exists(parsed.input_lookup):
		sys.exit('ERROR: %s does not exist.' % parsed.input_lookup)
	if not os.path.exists(parsed.gene_list):
		sys.exit('ERROR: %s does not exist.' % parsed.gene_list)
	
	count_matrix = create_count_matrix(parsed.input_lookup, parsed.gene_list)
	np.savetxt(parsed.output_count_file, count_matrix, delimiter=',', fmt='%s')
	

if __name__ == '__main__':
	main(sys.argv)
