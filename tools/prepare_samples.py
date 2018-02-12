#!/usr/bin/python
import sys
import argparse
import os.path
import pandas as pd
import glob
from utils import *


def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--metadata', required=True,
						help='Metadata file(s). Use delimiter "," if mutilple metadata files are used together for the same analysis group.')
	parser.add_argument('-s', '--samples', default='metadata/sample_summary.xlsx',
						help='Sample summary metadata file.')
	parser.add_argument('-g', '--group_num',
						help='Analysis group number. It is required if mutilple metadata files are used.')
	parser.add_argument('--condition_descriptors', default='TREATMENT,TIME_POINT',
						help='Experimental conditions to describe the sample. Use delimiter "," if multiple descriptors are used.')
	parser.add_argument('--qc_configure', default='tools/qc_config.yaml',
						help='Configuration file for quality assessment.')
	return parser.parse_args(argv[1:])


def build_dataframe(samples, metadata, df_cols, qc_cols, group):
	"""
	Build data frame for sample summary.
	"""
	df = initialize_dataframe(samples, df_cols, group)
	for md in metadata:
		if not os.path.exists(md):
			sys.exit('ERROR: %s does not exist.' % md)
		print '... Adding samples in %s' % md
		df = populate_dataframe(df, md, df_cols, qc_cols)
	## force group number if multiple metadata are given
	if group is not None:
		df['GROUP'] = group
	return df


def initialize_dataframe(samples, df_cols, group):
	"""
	Initialize the dataframe.
	"""
	if not os.path.exists(samples):
		print '... Creating sample summary %s' % samples
		df = pd.DataFrame(columns=df_cols)
	else:
		df = pd.read_excel(samples)
	return df


def populate_dataframe(df, metadata, df_cols, qc_cols, starting_sample=None):
	"""
	Add samples to dataframe from a metadata sheet.
	Find the corresponding fastq file.
	"""
	## read metadata
	df2 = pd.read_excel(metadata)
	df2 = df2.reset_index().drop(['index'], axis=1)
	# df2[qc_cols] = df2[qc_cols].apply(pd.to_numeric)
	if 'FILE' not in df2:
		df2 = pd.concat([df2, pd.Series([0]*df2.shape[0], name='FILE')], axis=1)
	sample = 0 if df.shape[0] == 0 else max(df['SAMPLE'])
	## update metadata
	for i,row in df2.iterrows():
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


def save_dataframe(filepath, df, df_cols):
	"""
	Save dataframe of quality assessment.
	"""
	df = df.sort_values(['SAMPLE'])
	if not filepath.endswith('.xlsx'):
		filepath += '.xlsx'
	df.to_excel(filepath, columns=df_cols, index=False)


def main(argv):
	parsed = parse_args(argv)
	## validate args
	metadata = [f.strip() for f in parsed.metadata.split(',')]
	if len(metadata) > 1 and parsed.group_num is None:
		sys.exit('ERROR: Group number is required.')

	## load QC config data
	global QC_dict
	QC_dict = load_config(parsed.qc_configure)

	## define columns in sample summary
	df_columns = ['GENOTYPE', 'SAMPLE', 'INDUCTION', 'LIBRARY', \
					'REPLICATE', 'INDEX', 'RUN_NUMBER', 'FILE', 'GROUP'] \
				+ [c.strip() for c in parsed.condition_descriptors.split(',')] \
				+ ['ST_PIPE', 'ST_TOTAL', 'ST_COMPLEXITY', 'ST_MUT_FOW', \
					'ST_RC_FOM', 'ST_COV_MED', 'AUTO_AUDIT', 'MANUAL_AUDIT', \
					'USER', 'NOTE']
	## define pre-alignment criteria to use
	qc_columns = ['BA', 'TAPESTATION', 'SPIKEIN_READS', 'SAMPLE_READS']
	## build sample summary sheet
	df = build_dataframe(parsed.samples, metadata, df_columns, qc_columns, parsed.group_num)
	save_dataframe(parsed.samples, df, df_columns)
	

if __name__ == '__main__':
	main(sys.argv)
