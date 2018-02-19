#!/usr/bin/python
import sys
import argparse
import os.path
import pandas as pd
from utils import *


def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-q', '--sample_quality', required=True,
						help='Sample quality assessment file.')
	parser.add_argument('-s', '--sample_summary', 
						default='metadata/sample_summary.xlsx',
						help='Sample summary metadata file.')
	parser.add_argument('--qc_configure', default='tools/qc_config.yaml',
						help='Configuration file for quality assessment.')
	return parser.parse_args(argv[1:])


def update_sample_summary(qa, samples):
	"""
	Decode status bit and update audit status afte manual audit of sample quality
	"""
	for i,row in qa.iterrows():
		sample_id = row['SAMPLE']
		## decompose status and update QA metric status
		decomp = decompose_status2bit(row['STATUS'])
		if decomp is not None:
			metrics = [QC_dict_rev[2**bit] for bit in decomp]
			for m in metrics:
				m_col = '_'.join(['ST', m])
				samples.loc[samples['SAMPLE'] == sample_id, m_col] = 1
		## update audit
		audit_cols = ['AUTO_AUDIT', 'MANUAL_AUDIT', 'USER', 'NOTE']
		for col in audit_cols:
			samples.loc[samples['SAMPLE'] == sample_id, col] = row[col]
	return samples


def save_dataframe(filepath, df, df_cols=None):
	"""
	Save dataframe of quality assessment.
	"""
	if not filepath.endswith('.xlsx'):
		filepath += '.xlsx'
	df.to_excel(filepath, index=False, columns=df_cols)


def reverse_QC_dict(QC_dict):
	"""
	A reversed dictionary with (status, metric) as (key, value)
	"""
	QC_dict_rev = {}
	for k in QC_dict.keys():
		if 'status' in QC_dict[k]:
			QC_dict_rev[QC_dict[k]['status']] = k
		else:
			v = QC_dict[k]
			if len(v.keys()) > 1:
				QC_dict_rev[v[v.keys()[0]]['status']] = k
	return QC_dict_rev


def main(argv):
	parsed = parse_args(argv)
	## validate args
	if not os.path.exists(parsed.sample_quality):
		sys.exit('ERROR: %s does not exist.' % parsed.sample_quality)
	if not os.path.exists(parsed.sample_summary):
		sys.exit('ERROR: %s does not exist.' % parsed.sample_summary)

	global QC_dict, QC_dict_rev
	QC_dict = load_config(parsed.qc_configure)
	QC_dict_rev = reverse_QC_dict(QC_dict)

	print '... Updating sample summary'
	qa = pd.read_excel(parsed.sample_quality)
	samples = pd.read_excel(parsed.sample_summary)
	df_columns = samples.columns.values
	samples = update_sample_summary(qa, samples)
	save_dataframe(parsed.sample_summary, samples, df_columns)


if __name__ == '__main__':
	main(sys.argv)