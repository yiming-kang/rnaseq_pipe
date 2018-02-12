#!/usr/bin/python
import sys
import argparse
import re
import pandas as pd
import numpy as np
from utils import *


def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-g', '--group_num', required=True, type=int, 
						help='Experiment group number.')
	parser.add_argument('-r', '--max_replicates', required=True, type=int,
						help='Maximal number of replicate in experiment design.')
	parser.add_argument('-w', '--wildtype', required=True,
						help='Wildtype genotype, e.g. CNAG_00000 for crypto, BY4741 for yeast.')
	parser.add_argument('-l', '--gene_list', required=True,
						help='Gene list.')
	parser.add_argument('-c', '--resistance_cassettes', required=True,
						help='Resistance cassettes inserted to replace the deleted genes. Use "," as delimiter if multiple cassettes exist.')
	parser.add_argument('-o', '--output_filepath', required=True,
						help='Filepath of sample quality summary.')
	parser.add_argument('-s', '--samples', default='metadata/sample_summary.xlsx',
						help='Sample summary metadata file.')
	parser.add_argument('--qc_configure', default='tools/qc_config.yaml',
						help='Configuration file for quality assessment.')
	parser.add_argument('--auto_audit_threshold', type=int, default=0,
						help='Threshold for automatical sample audit.')
	return parser.parse_args(argv[1:])


def initialize_dataframe(samples, df_cols, group):
	"""
	Define the QC dataframe.
	"""
	df = pd.DataFrame(columns=df_cols)
	df2 = pd.read_excel(samples, dtype=np.str)
	df2 = df2[df2['GROUP'] == group][['GENOTYPE','REPLICATE','SAMPLE']]
	df2 = df2.reset_index().drop(['index'], axis=1)
	df2 = pd.concat([df2, pd.Series([0]*df2.shape[0], name='STATUS')], axis=1)
	df2 = pd.concat([df2, pd.Series([0]*df2.shape[0], name='AUTO_AUDIT')], axis=1)
	return df.append(df2)


def combined_expression_data(df, gids, expr_tool='stringtie'):
	"""
	Combine individual expression profiles into a single expression matrix.
	Dimension = genes x samples.
	"""
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
	"""
	Assess percentage of uniquely mapped reads over all reads.
	"""
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
		row['TOTAL'] = total_reads
		row['COMPLEXITY'] = complexity
		if total_reads < QC_dict['TOTAL_READS']['threshold']:
			row['STATUS'] += QC_dict['TOTAL_READS']['status']
		if complexity < QC_dict['COMPLEXITY']['threshold']:
			row['STATUS'] += QC_dict['COMPLEXITY']['status']
		df.iloc[i] = row
	return df


def assess_efficient_mutation(df, expr, sample_dict, wt):
	"""
	Assess the completeness of gene deletion or efficiency of gene 
	overexpression by caluclating the expression of the perturbed gene
	in mutant sample over mean expression of the same gene in wildtype.
	"""
	## calculate mean expression level of each gene
	wt_expr = pd.Series(pd.DataFrame.mean(expr[sample_dict[wt].values()], 
						axis=1), name='mean_fpkm')
	wt_expr = pd.concat([expr, wt_expr], axis=1)
	## calculate efficiency of gene deletion, ignoring overexpression(*_over)
	for i,row in df[df['GENOTYPE'] != wt].iterrows():
		sample = row['GENOTYPE'] +'-'+ str(row['SAMPLE'])
		## check for each mutant gene (there could be multiple mutant genes, delimited by '.')
		mut_fow_list = []
		for mut_gene in row['GENOTYPE'].split('.'):
			mut_fow = float(expr[expr['gene'] == mut_gene][sample]) / \
					float(wt_expr[wt_expr['gene'] == mut_gene]['mean_fpkm'])
			if mut_gene.endswith('_over'):
				## check overexpression
				if (mut_fow < QC_dict['MUT_FOW']['OVEREXPRESSION']['threshold']) and (row['STATUS'] < QC_dict['MUT_FOW']['OVEREXPRESSION']['status']):
					row['STATUS'] += QC_dict['MUT_FOW']['OVEREXPRESSION']['status']
			else:
				## check deletion
				if (mut_fow > QC_dict['MUT_FOW']['DELETION']['threshold']) and (row['STATUS'] < QC_dict['MUT_FOW']['DELETION']['status']):
					row['STATUS'] += QC_dict['MUT_FOW']['DELETION']['status']
			mut_fow_list.append(str(mut_fow))
		row['MUT_FOW'] = ','.join(mut_fow_list)
		df.iloc[i] = row
	return df


def assess_resistance_cassettes(df, expr, resi_cass, wt):
	"""
	Assess drug resistance marker gene expression, making sure the proper
	marker gene is swapped in place of the perturbed gene.
	"""
	## get the median of resistance cassettes
	rc_med_dict = {}
	mut_samples = [s for s in expr.columns.values if (not s.startswith(wt)) and (s != 'gene')]
	for rc in resi_cass:
		## exclude wildtypes and samples that have other makers expressed 
		rc0 = [x for x in resi_cass if x != rc]
		rc0_fpkm = expr.loc[expr['gene'].isin(rc0), mut_samples]
		rc_fpkm = expr.loc[expr['gene'] == rc, mut_samples]
		rc_fpkm = rc_fpkm.loc[:, (np.sum(rc_fpkm, axis=0) != 0) & (np.sum(rc0_fpkm, axis=0) == 0)]
		rc_med_dict[rc] = rc_fom = np.nan if rc_fpkm.empty else np.median(rc_fpkm)
	## calcualte FOM (fold change over mutant) of the resistance cassette
	for i,row in df.iterrows():
		genotype = row['GENOTYPE']
		sample = genotype +'-'+ str(row['SAMPLE'])
		## update FOM
		for rc in rc_med_dict.keys():
			row[rc+'_FOM'] = np.nan if np.isnan(rc_med_dict[rc]) else float(expr.loc[expr['gene'] == rc, sample])/rc_med_dict[rc]
		## flag those two problems: 
		## 1. the resistance cassette is exprssed in WT
		## 2. more than one resistance cassette is expressed in single mutant
		## TODO: add criteria for multi-mutants 
		fom_check = [row[rc+'_FOM'] > QC_dict['MARKER_FOM']['threshold']*rc_med_dict[rc] for rc in rc_med_dict.keys()]
		if genotype == wt and sum(fom_check) > 0:
			row['STATUS'] += QC_dict['MARKER_FOM']['status']
		if genotype != wt and len(genotype.split('.')) > 1 and sum(fom_check) > 1:
			row['STATUS'] += QC_dict['MARKER_FOM']['status']
	return df


def assess_replicate_concordance(df, expr, sample_dict):
	"""
	Assess the concordance among the replicates of each genotype by calculating
	the COV of each combination of replicates. Then find the maximal number of
	concordant replicates.
	"""
	cov = expr['gene']
	## calcualte COV medians for replicate combinations
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
		## find the maximal number of replicates that pass concordance threshold
		for rep_num in sorted(cov_meds_dict.keys())[::-1]:
			rep_combo = cov_meds_dict[rep_num]['rep_combos']
			cov_meds = cov_meds_dict[rep_num]['cov_meds']
			if sum([c < QC_dict['COV_MED']['threshold'] for c in cov_meds]) > 0:
				best_combo = rep_combo[np.argmin(cov_meds)]
				break 
		max_rep_combo = cov_meds_dict[max(cov_meds_dict.keys())]['rep_combos'][0]
		outlier_reps = set(max_rep_combo) - set(best_combo)
		for rep in outlier_reps:
			df.loc[(df['GENOTYPE'] == genotype) & (df['REPLICATE'] == rep), 'STATUS'] += QC_dict['COV_MED']['status']
	return df


def calculate_cov_median(x):
	"""
	Calculate the median of COVs (coefficient of variation) among replicates
	"""
	covs = np.std(x, axis=1) / np.mean(x, axis=1)
	return np.nanmedian(covs)


def update_auto_audit(df, threshold):
	"""
	Automatically flag sample with status over threshold 
	"""
	df['AUTO_AUDIT'][np.where(df['STATUS'] > threshold)[0]] = 1
	return df


def save_dataframe(filepath, df, df_cols):
	"""
	Save dataframe of quality assessment
	"""
	df = df.sort_values(['GENOTYPE','REPLICATE'])
	if not filepath.endswith('.xlsx'):
		filepath += '.xlsx'
	df.to_excel(filepath, columns=df_cols, index=False, freeze_panes=(1,3))


def main(argv):
	parsed = parse_args(argv)
	## validate args
	if os.path.exists(parsed.output_filepath):
		sys.exit('WARNING: %s already exists, rename the file to proceed.' % parsed.output_filepath)
	if not os.path.exists(parsed.samples):
		sys.exit('ERROR: %s does not exist.' % parsed.samples)
	if not os.path.exists(os.path.dirname(parsed.output_filepath)):
		sys.exit('ERROR: %s does not exist.' % os.path.dirname(parsed.output_filepath))

	## load QC config data
	## TODO: complexity.thresh <- mean(alignment.sum$COMPLEXITY[indx]) - 2*sd(alignment.sum$COMPLEXITY[indx]);
	global QC_dict
	QC_dict = load_config(parsed.qc_configure)

	## do QA
	print '... Preparing QA dataframe'
	resistance_cassettes = [rc.strip() for rc in parsed.resistance_cassettes.split(',')]
	df_columns = ['GENOTYPE','REPLICATE','SAMPLE', \
					'TOTAL','COMPLEXITY','MUT_FOW'] \
				+ [rc+'_FOM' for rc in resistance_cassettes] \
				+ ['COV_MED_REP'+''.join(np.array(combo, dtype=str)) for combo in make_combinations(range(1,parsed.max_replicates+1))] \
				+ ['STATUS', 'AUTO_AUDIT', 'MANUAL_AUDIT', 'USER', 'NOTE']
	df = initialize_dataframe(parsed.samples, df_columns, parsed.group_num)
	expr, sample_dict = combined_expression_data(df, parsed.gene_list)
	print '... Assessing reads mapping'
	df = assess_mapping_quality(df)
	print '... Assessing efficiency of gene mutation'
	df = assess_efficient_mutation(df, expr, sample_dict, parsed.wildtype)
	print '... Assessing insertion of resistance cassette'
	df = assess_resistance_cassettes(df, expr, resistance_cassettes, parsed.wildtype)
	print '... Assessing concordance among replicates'
	df = assess_replicate_concordance(df, expr, sample_dict)
	df = update_auto_audit(df, parsed.auto_audit_threshold)
	save_dataframe(parsed.output_filepath, df, df_columns)
	

if __name__ == '__main__':
	main(sys.argv)
