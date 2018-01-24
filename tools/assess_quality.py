#!/usr/bin/python
import sys
import argparse
import re
import pandas as pd


global QC_dict
QC_dict = {'total_reads': {'threshold': 5*(10**6), 'status': 1},
			'complexity': {'threshold': .1, 'status': 2},
			'delfow': {'threshold': .1, 'status': 4}}


def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--samples',
						help='Sample summary metadata file.')
	parser.add_argument('-e', '--experiment_group', type=int,
						help='Identifier for the experiment group.')
	parser.add_argument('-w', '--wildtype',
						help='Wildtype genotyp, e.g. CNAG_00000 for crypto, BY4741 for yeast.')
	parser.add_argument('-g', '--gene_list',
						help='Gene list.')
	return parser.parse_args(argv[1:])


def initialize_dataframe(samples, group):
	## define dataframe
	df = pd.DataFrame(columns=['GENOTYPE','REPLICATE','SAMPLE','TOTAL_READS','COMPLEXITY','DELFOW','STATUS'])
	## load sample metadata
	df2 = pd.read_csv(samples, delimiter='\t')
	df2 = df2[df2['GROUP'] == group][['GENOTYPE','REPLICATE','SAMPLE']]
	## initialize status
	df2 = pd.concat([df2, pd.Series([0]*df2.shape[0], name='STATUS')], axis=1)
	return df.append(df2)


def combined_expression_data(df, gids, expr_tool='stringtie'):
	## combine into a expression matrix, genes x samples
	expr = pd.read_csv(gids, names=['gene'])
	for i,row in df.iterrows():
		## get expression profile
		sample = row['GENOTYPE'] +'-'+ str(row['SAMPLE'])
		filepath = '/'.join(['expression', expr_tool+'_fpkm', sample+'.expr'])
		indiv_expr = pd.read_csv(filepath, names=[sample])
		## concatenate horizontally
		expr = pd.concat([expr, indiv_expr], axis=1)
	return expr


def assess_mapping_quality(df, aligner_tool='novoalign'):
	for i,row in df.iterrows():
		sample = row['GENOTYPE'] +'-'+ str(row['SAMPLE'])
		filepath = '/'.join(['alignment', aligner_tool, sample, aligner_tool+'.log'])
		## read alignment log
		reader = open(filepath, 'r')
		lines = reader.readlines()
		for line in lines:
			reg_total = re.search(r'Read Sequences: (.\d+)', line)
			reg_uniq = re.search(r'Unique Alignment: (.\d+)', line)
			if reg_total:
				total_reads = int(reg_total.group(1))
			if reg_uniq:
				uniq_mapped_reads = int(reg_uniq.group(1))
		reader.close()
		complexity = uniq_mapped_reads/float(total_reads)
		## set mapping quality
		row['TOTAL_READS'] = total_reads
		row['COMPLEXITY'] = complexity
		if total_reads < QC_dict['total_reads']['threshold']:
			row['STATUS'] += QC_dict['total_reads']['status']
		if complexity < QC_dict['complexity']['threshold']:
			row['STATUS'] += QC_dict['complexity']['status']
		df.iloc[i] = row
	return df


def assess_efficient_perturbation(df, expr, wt):
	## get gene levels in wildtype samples
	wt_samples = []
	for i,row in df[df['GENOTYPE'] == wt].iterrows():
		sample = row['GENOTYPE'] +'-'+ str(row['SAMPLE'])
		wt_samples.append(sample)
	## calculate mean expression level of each gene
	wt_expr = pd.Series(pd.DataFrame.mean(expr[wt_samples], 
					axis=1), name='mean_fpkm')
	wt_expr = pd.concat([expr, wt_expr], axis=1)
	## calculate efficiency of marker gene deletion, ignoring overexpression(*_over)
	for i,row in df[df['GENOTYPE'] != wt].iterrows():
		sample = row['GENOTYPE'] +'-'+ str(row['SAMPLE'])
		## check for each marker gene (there could be multiple marker gene perturbation, delimited by '.')
		ineff_pertub = False
		delfow_list = []
		for marker_gene in row['GENOTYPE'].split('.'):
			if not marker_gene.endswith('_over'):
				marker_gene_pert = float(expr[expr['gene'] == marker_gene][sample])
				marker_gene_wt = float(wt_expr[wt_expr['gene'] == marker_gene]['mean_fpkm'])
				delfow = marker_gene_pert/marker_gene_wt
				if delfow > QC_dict['delfow']['threshold']:
					ineff_pertub = True
				delfow_list.append(str(delfow))
		## set inefficient deletion
		row['DELFOW'] = ','.join(delfow_list)
		if ineff_pertub:
			row['STATUS'] += QC_dict['delfow']['status']
		df.iloc[i] = row
	return df


def assess_relicate_concordance(df, expr):
	pass


def main(argv):
	parsed = parse_args(argv)

	df = initialize_dataframe(parsed.samples, parsed.experiment_group)
	expr = combined_expression_data(df, parsed.gene_list)
	df = assess_mapping_quality(df)
	df = assess_efficient_perturbation(df, expr, parsed.wildtype)
	df = assess_replicate_concordance(df, expr)
	print df

if __name__ == '__main__':
	main(sys.argv)
