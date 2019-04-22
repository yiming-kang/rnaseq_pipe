#!/usr/bin/python
import sys
import os.path
import pandas as pd
import argparse

def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--sample_summary', required=True, 
						help='Sample metadata sheet after QA auditing.')
	parser.add_argument('-o', '--report_file', required=True,
						help='Sample summary report.')
	return parser.parse_args(argv[1:])


def analyze_summary(df):
	report = pd.DataFrame({"GENOTYPE":[], "TOTAL":[], "PASS":[], "FAIL":[], "NOTE":[], "ID":[]})
	## iterate thru to curate stats
	for _,row in df.iterrows():
		g = row["GENOTYPE"] 
		i = row["SAMPLE"]
		a = row["MANUAL_AUDIT"]
		n = row["NOTE"]
		p = row["ST_PIPE"]
		if g not in report["GENOTYPE"].values:
			report = report.append(pd.DataFrame({"GENOTYPE":[g], "TOTAL":[0], "PASS":[0], "FAIL":[0], "NOTE":[""], "ID":[""]}))
		## count total reps
		report.loc[report["GENOTYPE"]==g, "TOTAL"] += 1
		## curate sample IDs
		report.loc[report["GENOTYPE"]==g, "ID"] += " | " + str(i)
		if a == 0:
			## curate pass samples
			report.loc[report["GENOTYPE"]==g, "PASS"] += 1
		else:
			## curate user's note on the failure reasons
			report.loc[report["GENOTYPE"]==g, "FAIL"] += 1
			report.loc[report["GENOTYPE"]==g, "NOTE"] += " | "+str(n)
	return report.sort_values(by=["PASS","GENOTYPE"])


def main(argv):
	parsed = parse_args(argv)
	## validate args
	if not os.path.exists(parsed.sample_summary):
		sys.exit('ERROR: %s does not exist.' % parsed.sample_summary)

	## load sample summary
	sample_summary = pd.read_excel(parsed.sample_summary)
	## analyze summary
	report = analyze_summary(sample_summary)
	## write report
	report.to_excel(parsed.report_file, columns=["GENOTYPE","TOTAL","PASS","FAIL","NOTE","ID"], index=False)


if __name__ == '__main__':
	main(sys.argv)
