#!/usr/bin/python
import sys
import argparse
import os.path
import pandas as pd
import glob
import copy
from utils import *


def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--samples', required=True,
						help='Sample summary metadata file.')
	parser.add_argument('-g', '--group_num', required=True,
						help='Analysis group number. It is required if mutilple metadata files are used.')
	parser.add_argument('-d', '--design_table', required=True,
						help='Auto-gen design table for DE anlaysis.')
	parser.add_argument('-w',  '--wildtype',
						help='Wildtype genotype, e.g. CNAG_00000 for crypto, BY4741 for yeast.')
	parser.add_argument('--condition_descriptors', default='TREATMENT,TIME_POINT',
						help='Experimental conditions to describe the sample. Use delimiter "," if multiple descriptors are used. Default is TREATMENT,TIME_POINT')
	return parser.parse_args(argv[1:])


def build_design_table(summary_df, cmp_cols, wt=None):
	"""
	Automatically build design table that has basic comparison groups.
	"""
	print '... Building design table'
	## intialize design dataframe
	design_df = summary_df[['GENOTYPE', 'REPLICATE', 'SAMPLE'] + cmp_cols].copy()
	cmp_cols = ['GENOTYPE'] + cmp_cols
	## create a dict for col: vals
	col_dict = {}
	for col in cmp_cols:
		col_dict[col] = pd.unique(summary_df[col])
	## iterate thru columns to compare
	for col in cmp_cols:
		if len(col_dict[col]) < 2: 
			## nothing to contrast with
			continue	
		## get value combination of other column descriptors
		other_cols = copy.copy(cmp_cols)
		other_cols.remove(col)
		n_other_cols = len(other_cols)
		other_vals = [col_dict[ocol] for ocol in other_cols]
		other_vals_combos = make_list_product(other_vals)
		## iterate thru with other column value fixed
		for vcombo in other_vals_combos:
			## get rows matching the column value combo
			df2 = summary_df[sum([summary_df[other_cols[k]] == vcombo[k] \
								for k in range(n_other_cols)]) == n_other_cols]
			vcombo_name = '-'.join([':'.join([other_cols[k], vcombo[k]]) \
									for k in range(n_other_cols)])
			## assess diff column types
			if col == 'GENOTYPE' and wt is not None:
				vals = col_dict[col].tolist()
				vals.remove(wt)
				for val in vals:
					## concatenate new column and set flag
					new_col = ''.join(['[',vcombo_name,']',col,':',wt,'-',val])
					design_df = pd.concat([design_df, pd.Series(['']*design_df.shape[0], name=new_col)], axis=1)
					design_df.loc[df2.index[df2[col] == wt], new_col] = '0'
					design_df.loc[df2.index[df2[col] == val], new_col] = '1'
			else: 
				## TODO: extend the function for other columns
				pass
	## sort dataframe
	design_df = design_df.sort_values(cmp_cols + ['SAMPLE'])
	return design_df


def save_dataframe(filepath, df, df_cols=None):
	"""
	Save dataframe of quality assessment.
	"""
	if not filepath.endswith('.xlsx'):
		filepath += '.xlsx'
	df.to_excel(filepath, index=False, columns=df_cols)



def main(argv):
	parsed = parse_args(argv)
	if not os.path.exists(parsed.samples):
		sys.exit('ERROR: %s does not exist.' % parsed.samples)

	## get conditions
	conditions = [c.strip() for c in parsed.condition_descriptors.split(',')]

	## load sample summary
	summary_df = pd.read_excel(parsed.samples)
	## prepare design table
	design_df = build_design_table(summary_df[summary_df['GROUP']==parsed.group_num], conditions, parsed.wildtype)
	save_dataframe(parsed.design_table, design_df)


if __name__ == '__main__':
	main(sys.argv)