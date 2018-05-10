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
	parser.add_argument('-e', '--expr_matrix', required=True,
						help='Input expression matrix.')
	parser.add_argument('--condition_descriptors', default='TREATMENT,TIME_POINT',
						help='Experimental conditions that describe the sample are used to identify subgroups within each genotype. Use delimiter "," if multiple descriptors are used.')
	parser.add_argument('-o', '--output_matrix',
						help='Output count matrix in csv format.')
	return parser.parse_args(argv[1:])


def query_metadata_features(metadata, descriptors):
	## query user input
	features = ['GENOTYPE'] + descriptors + ['INDUCTION', 'LIBRARY']
	query_dict = {}
	for feature in features:
		## print messages for querying user input
		feature_values = [x.encode("utf-8") for x in pd.unique(metadata[feature])]
		print "... Querying %s :" % feature  
		print np.array(feature_values)
		user_input = raw_input("Select from above: (Press 'Enter' to select all. If mutliple values selected, use delimiter ',')\n")
		## parse user input
		if user_input.strip() == "":
			query_dict[feature] = feature_values
		else:
			query_dict[feature] = [x.strip() for x in user_input.split(",")]
	## filter metadata based on query
	sub_meta = metadata.copy()
	for feature, values in query_dict.iteritems():
		sub_meta = sub_meta.loc[sub_meta[feature].isin(values),]
	sub_meta = sub_meta.sort_values(features+['SAMPLE'])
	## get sample list
	queried_samples = []
	for i,row in sub_meta[['GENOTYPE','SAMPLE']].iterrows():
		row = [str(x).encode("utf-8") for x in row]
		queried_samples.append('-'.join(row))
	return queried_samples


def query_genes_of_interest():
	## print messages for querying user input
	print "... Querying genes of interest :"
	user_input = raw_input("Type in genes of interest: (Press 'Enter' to select all. If mutliple values selected, use delimiter ',')\n")
	## parse user input
	if user_input.strip() == "":
		queried_genes = None
	else:
		queried_genes = [x.strip() for x in user_input.split(",")]
	return queried_genes


def main(argv):
	parsed = parse_args(argv)
	if not os.path.exists(parsed.samples):
		sys.exit('ERROR: %s does not exist.' % parsed.samples)
	if not os.path.exists(parsed.expr_matrix):
		sys.exit('ERROR: %s does not exist.' % parsed.expr_matrix)
	condition_descriptors = [x.strip() for x in parsed.condition_descriptors.split(',')]
	## query samples and genes
	metadata = pd.read_excel(parsed.samples)
	sample_list = query_metadata_features(metadata, condition_descriptors)
	gene_list = query_genes_of_interest()
	## load expression matrix and subset based on the query
	expr_matrix = pd.read_csv(parsed.expr_matrix, index_col=0)
	output_matrix = expr_matrix.loc[gene_list, np.intersect1d(sample_list, expr_matrix.columns.values)]
	print "\n... Expression matrix filtered from query:"
	print output_matrix
	## optionally save output matrix
	if parsed.output_matrix:
		np.savetxt(parsed.output_matrix, output_matrix, delimiter=',', fmt='%s')
	

if __name__ == '__main__':
	main(sys.argv)