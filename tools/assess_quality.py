#!/usr/bin/python
import sys
import argparse
import re
import pandas as pd
import numpy as np
from itertools import combinations


global QC_dict
QC_dict = {'total_reads': {'threshold': 5*(10**6), 'status': 1},
			'complexity': {'threshold': .1, 'status': 2},
			'delfow': {'threshold': .1, 'status': 4},
			'cov_med': {'threshold': .25, 'status': 8}}


def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--samples',
						help='Sample summary metadata file.')
	parser.add_argument('-e', '--experiment_group', type=int,
						help='Identifier for the experiment group.')
	parser.add_argument('-w', '--wildtype',
						help='Wildtype genotype, e.g. CNAG_00000 for crypto, BY4741 for yeast.')
	parser.add_argument('-g', '--gene_list',
						help='Gene list.')
	parser.add_argument('-o', '--output_filepath',
						help='Filepath of sample quality summary.')
	parser.add_argument('-r', '--max_replicates', default=4,
						help='Maximal number of replicate in experiment design.')
	return parser.parse_args(argv[1:])


def initialize_dataframe(samples, df_cols, group):
	## define dataframe
	df = pd.DataFrame(columns=df_cols)
	df2 = pd.read_csv(samples, delimiter='\t')
	df2 = df2[df2['GROUP'] == group][['GENOTYPE','REPLICATE','SAMPLE']]
	## initialize status
	df2 = pd.concat([df2, pd.Series([0]*df2.shape[0], name='STATUS')], axis=1)
	return df.append(df2)


def combined_expression_data(df, gids, expr_tool='stringtie'):
	## combine into a expression matrix, genes x samples
	expr = pd.read_csv(gids, names=['gene'])
	sample_dict = {}
	for i,row in df.iterrows():
		## get expression profile
		genotype = row['GENOTYPE']
		sample = genotype +'-'+ str(row['SAMPLE'])
		filepath = '/'.join(['expression', expr_tool+'_fpkm', sample+'.expr'])
		indiv_expr = pd.read_csv(filepath, names=[sample])
		## concatenate horizontally
		expr = pd.concat([expr, indiv_expr], axis=1)
		## update sample dictionary
		if genotype not in sample_dict.keys():
			sample_dict[genotype] = {}
		sample_dict[genotype][row['REPLICATE']] = sample
	return expr, sample_dict


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


def assess_efficient_perturbation(df, expr, sample_dict, wt):
	## calculate mean expression level of each gene
	wt_expr = pd.Series(pd.DataFrame.mean(expr[sample_dict[wt].values()], 
					axis=1), name='mean_fpkm')
	wt_expr = pd.concat([expr, wt_expr], axis=1)
	## calculate efficiency of marker gene deletion, ignoring overexpression(*_over)
	for i,row in df[df['GENOTYPE'] != wt].iterrows():
		sample = row['GENOTYPE'] +'-'+ str(row['SAMPLE'])
		## check for each marker gene (there could be multiple marker gene perturbation, delimited by '.')
		delfow_list = []
		for marker_gene in row['GENOTYPE'].split('.'):
			if not marker_gene.endswith('_over'):
				marker_gene_pert = float(expr[expr['gene'] == marker_gene][sample])
				marker_gene_wt = float(wt_expr[wt_expr['gene'] == marker_gene]['mean_fpkm'])
				delfow = marker_gene_pert/marker_gene_wt
				## set inefficient deletion
				if (delfow > QC_dict['delfow']['threshold']) and (row['STATUS'] < QC_dict['delfow']['status']):
					row['STATUS'] += QC_dict['delfow']['status']
				delfow_list.append(str(delfow))
		row['DELFOW'] = ','.join(delfow_list)
		df.iloc[i] = row
	return df


def assess_replicate_concordance(df, expr, sample_dict):
	cov = expr['gene']
	## calcualte COV medians for each genotype's replicate combinations
	for genotype in sorted(sample_dict.keys()):
		cov_meds_dict = {}
		rep_combos = make_combinations(sample_dict[genotype].keys())
		for rep_combo in rep_combos:
			sample_combo = [sample_dict[genotype][rep] for rep in sorted(rep_combo)]
			## calculate COV median
			cov_median = calculate_cov_median(expr[sample_combo])
			rep_combo_col = 'COV_MED_REP'+''.join(np.array(rep_combo, dtype=str))
			df.loc[df['GENOTYPE'] == genotype, rep_combo_col] = cov_median
			## store COV median at the respective rep number
			rep_num = len(rep_combo)
			if rep_num not in cov_meds_dict.keys():
				cov_meds_dict[rep_num] = {'rep_combos': [], 'cov_meds': []}
			cov_meds_dict[rep_num]['rep_combos'].append(rep_combo)
			cov_meds_dict[rep_num]['cov_meds'].append(cov_median)
		## assess replicate concordance
		for rep_num in sorted(cov_meds_dict.keys())[::-1]:
			rep_combo = cov_meds_dict[rep_num]['rep_combos']
			cov_meds = cov_meds_dict[rep_num]['cov_meds']
			if sum([c < QC_dict['cov_med']['threshold'] for c in cov_meds]) > 0:
				best_combo = rep_combo[np.argmin(cov_meds)]
				break 
		max_rep_combo = cov_meds_dict[max(cov_meds_dict.keys())]['rep_combos'][0]
		outlier_reps = set(max_rep_combo) - set(best_combo)
		for rep in outlier_reps:
			df.loc[(df['GENOTYPE'] == genotype) & (df['REPLICATE'] == rep), 'STATUS'] += QC_dict['cov_med']['status']
	return df


def make_combinations(lst):
	## make all possible replicate combinations
	if len(lst) < 2:
		return [lst]
	combo = []
	for i in range(len(lst), 1, -1):
		for s in combinations(lst, i):
			combo.append(s)
	return combo


def calculate_cov_median(x):
	## 
	covs = np.std(x, axis=1) / np.mean(x, axis=1)
	return np.nanmedian(covs)


def save_dataframe(filepath, df, df_cols):
	df = df.sort_values(['GENOTYPE','REPLICATE'])
	df.to_csv(filepath, sep='\t', columns=df_cols, index=False, float_format='%.3f')


def main(argv):
	parsed = parse_args(argv)

	df_columns = ['GENOTYPE','REPLICATE','SAMPLE',
				'TOTAL_READS','COMPLEXITY','DELFOW'] + \
				['COV_MED_REP'+''.join(np.array(combo, dtype=str)) \
				for combo in make_combinations(range(1,parsed.max_replicates+1))] + \
				['STATUS']
	df = initialize_dataframe(parsed.samples, df_columns, parsed.experiment_group)
	expr, sample_dict = combined_expression_data(df, parsed.gene_list)
	df = assess_mapping_quality(df)
	df = assess_efficient_perturbation(df, expr, sample_dict, parsed.wildtype)
	df = assess_replicate_concordance(df, expr, sample_dict)
	save_dataframe(parsed.output_filepath, df, df_columns)

if __name__ == '__main__':
	main(sys.argv)
