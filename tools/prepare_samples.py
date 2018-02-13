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
	parser.add_argument('-m', '--metadata', required=True,
						help='Metadata file(s). Use delimiter "," if mutilple metadata files are used together for the same analysis group.')
	parser.add_argument('-g', '--group_num', required=True, type=int,
						help='Analysis group number. It is required if mutilple metadata files are used.')
	parser.add_argument('-s', '--samples', default='metadata/sample_summary.xlsx',
						help='Sample summary metadata file.')
	parser.add_argument('-d', '--design_table',
						help='Auto-gen design table for DE anlaysis.')
	parser.add_argument('-w',  '--wildtype',
						help='Wildtype genotype, e.g. CNAG_00000 for crypto, BY4741 for yeast.')
	parser.add_argument('--condition_descriptors', default='TREATMENT,TIME_POINT',
						help='Experimental conditions to describe the sample. Use delimiter "," if multiple descriptors are used.')
	parser.add_argument('--qc_configure', default='tools/qc_config.yaml',
						help='Configuration file for quality assessment.')
	return parser.parse_args(argv[1:])


def build_sample_summary(samples, metadata, df_cols, qc_cols, group, conditions):
	"""
	Build data frame for sample summary.
	"""
	df = initialize_sample_summary(samples, df_cols)
	for md in metadata:
		if not os.path.exists(md):
			sys.exit('ERROR: %s does not exist.' % md)
		print '... Adding samples in %s' % md
		df = populate_sample_summary(df, md, df_cols, qc_cols, conditions)
	## force group number if multiple metadata are given
	if group is not None:
		df['GROUP'] = group
	return df


def initialize_sample_summary(samples, df_cols):
	"""
	Initialize the dataframe.
	"""
	if not os.path.exists(samples):
		print '... Creating sample summary %s' % samples
		df = pd.DataFrame(columns=df_cols)
	else:
		df = pd.read_excel(samples)
	return df


def populate_sample_summary(df, metadata, df_cols, qc_cols, conditions):
	"""
	Add samples to dataframe from a metadata sheet.
	Find the corresponding fastq file.
	"""
	## get exisiting samples, used to check adding redudant samples
	check_cols = ['GENOTYPE', 'REPLICATE', 'INDUCTION', 'LIBRARY'] + conditions
	exist_samples = get_exisiting_samples(df, check_cols)
	## read metadata
	df2 = pd.read_excel(metadata)
	df2 = df2.reset_index().drop(['index'], axis=1)
	# df2[qc_cols] = df2[qc_cols].apply(pd.to_numeric)
	if 'FILE' not in df2:
		df2 = pd.concat([df2, pd.Series([0]*df2.shape[0], name='FILE')], axis=1)
	sample = 0 if df.shape[0] == 0 else max(df['SAMPLE'])
	## update metadata
	for i,row in df2.iterrows():
		## check sample redudancy
		sample_descriptor = [row[col] for col in sorted(check_cols)]
		if sample_descriptor in exist_samples:
			sys.exit('WARNING: found existing sample:\n%s\nCheck the redudancy issue.' % ' '.join([':'.join([col, str(row[col])]) for col in sorted(check_cols)]))
		## add the unique sample id
		sample += 1
		df2.loc[i, 'SAMPLE'] = sample
		## add audit flag on pre-alignment QC
		if sum([row[qc] < QC_dict[qc]['threshold'] for qc in qc_cols]) != 0:
			df2.loc[i, 'ST_PIPE'] = 1
		## valid uniqueness and existence of the fastq file
		file = 'sequence/run_'+ str(row['RUN_NUMBER']) +'_samples/*'+ row['INDEX'] +'*'
		files_found = glob.glob(file)
		if len(files_found) == 0:
			sys.exit('ERROR: No file exists with index %s and run number %s.' % (row['INDEX'], row['RUN_NUMBER']))
		elif len(files_found) > 1:
			sys.exit('ERROR: Mutilple files exist with the same index %s and run number %s.' % (row['INDEX'], row['RUN_NUMBER']))
		else:
			df2.loc[i, 'FILE'] = files_found[0]
	df = df.append(df2, ignore_index=True)[df_cols]
	return df


def get_exisiting_samples(df, check_cols):
	"""
	Get exisiting samples with respective descriptors.
	"""
	samples = []
	for i,row in df.iterrows():
		sample_descriptor = [row[col] for col in sorted(check_cols)]
		samples.append(sample_descriptor)
	return samples


def save_dataframe(filepath, df, df_cols=None):
	"""
	Save dataframe of quality assessment.
	"""
	df = df.sort_values(['SAMPLE'])
	if not filepath.endswith('.xlsx'):
		filepath += '.xlsx'
	if df_cols is None:
		df.to_excel(filepath, index=False)
	else:
		df.to_excel(filepath, columns=df_cols, index=False)


def build_design_table(summary_df, cmp_cols, wt=None):
	"""
	Automatically build design table that has basic comparison groups.
	"""
	## intialize design dataframe
	design_df = summary_df[['GENOTYPE', 'REPLICATE']].copy()
	## create a dict for col: vals
	col_dict = {}
	for col in cmp_cols:
		col_dict[col] = pd.unique(summary_df[col])
	## iterate thru columns to compare
	for col in cmp_cols:
		if len(col_dict[col]) > 1:
			## get value combination of other column descriptors
			other_cols = copy.copy(cmp_cols)
			other_cols.remove(col)
			other_vals = [col_dict[ocol] for ocol in other_cols]
			other_vals_combo = make_list_product(other_vals)
			## assess diff column types
			if col == 'GENOTYPE' and wt is not None:
				col ## TODO
			else: ## TODO: extend the function for other columns
				pass
	return design_df


def main(argv):
	parsed = parse_args(argv)
	## validate args
	metadata = [f.strip() for f in parsed.metadata.split(',')]
	if len(metadata) > 1 and parsed.group_num is None:
		sys.exit('ERROR: Group number is required.')

	## load QC config data
	global QC_dict
	QC_dict = load_config(parsed.qc_configure)
	## get conditions
	conditions = [c.strip() for c in parsed.condition_descriptors.split(',')]

	## define columns in sample summary
	df_columns = ['GENOTYPE', 'SAMPLE', 'INDUCTION', 'LIBRARY', \
					'REPLICATE', 'INDEX', 'RUN_NUMBER', 'FILE', 'GROUP'] \
				+ conditions + ['ST_PIPE', 'ST_TOTAL', 'ST_COMPLEXITY', \
					'ST_MUT_FOW', 'ST_RC_FOM', 'ST_COV_MED', 'AUTO_AUDIT', \
					'MANUAL_AUDIT', 'USER', 'NOTE']
	## define pre-alignment criteria to use
	qc_columns = ['BA', 'TAPESTATION', 'SPIKEIN_READS', 'SAMPLE_READS']
	## build sample summary sheet
	summary_df = build_sample_summary(parsed.samples, metadata, df_columns, qc_columns, parsed.group_num, conditions)
	# save_dataframe(parsed.samples, summary_df, df_columns)

	## build design table
	comparison_columns = ['GENOTYPE'] + conditions
	design_df = build_design_table(df[df['GROUP']==parsed.group_num], comparison_columns, parsed.wildtype)
	# save_dataframe(parsed.design_table, design_df)


if __name__ == '__main__':
	main(sys.argv)
