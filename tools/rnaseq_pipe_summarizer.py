import sys
import os
import re
import argparse
from glob import glob
import pandas as pd


ALIGN_VARS = ["Read Sequences", "Unique Alignment", "Multi Mapped", "No Mapping Found"]
COUNT_VARS = ["total_mapped_reads", "with_feature", "no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique"]


def parse_args(argv):
	parser = argparse.ArgumentParser(description="This script summarizes the output from pipeline wrapper.")
	parser.add_argument("--align_path",
						help="Directory for alginment log files. This currently only works for Novoalign output.")
	parser.add_argument("--count_path",
						help="Directory for gene read count files. This currently only works for HT-Seq output.")
	parser.add_argument("--output_csv",
						help="File path for the pipeline summary.")
	args = parser.parse_args(argv[1:])
	return args


def compile_data(dir_path, suffix):
	df = pd.DataFrame()
	file_paths = glob("{}/*{}".format(dir_path, suffix))
	for file_path in file_paths:
		sample = re.findall(r'(.+?){0}'.format(suffix), os.path.basename(file_path))[0]
		if "novoalign" in suffix:
			data = parse_alignment_log(file_path)
		elif "read_count" in suffix:
			data = parse_gene_count(file_path)
		data.update({"Sample": sample})
		df = df.append(pd.Series(data), ignore_index=True)
	return df.set_index("Sample")


def parse_alignment_log(file_path):
	with open(file_path, "r") as f:
		lines = f.readlines()
	out = {}
	for line in lines:
		line = line.strip("#").strip()
		for k in ALIGN_VARS:
			if line.startswith(k):
				v = int(line.split(":")[-1].split("(")[0].strip())
				out[k] = v if k == ALIGN_VARS[0] else v/float(out[ALIGN_VARS[0]])
	return out


def parse_gene_count(file_path):
	with open(file_path, "r") as f:
		lines = f.readlines()
	out = {COUNT_VARS[1]: 0}
	for line in lines:
		k, v = line.strip().split()
		if line.startswith("__"):
			out[k[2:]] = int(v)
		else:
			out[COUNT_VARS[1]] += int(v)
	out[COUNT_VARS[0]] = sum([out[k] for k in COUNT_VARS[1:]])
	for k in COUNT_VARS[1:]:
		out[k] /= float(out[COUNT_VARS[0]])
	return out


def main(argv):
	args = parse_args(argv)
	align_path = args.align_path
	count_path = args.count_path
	output_csv = args.output_csv

	align_df = compile_data(align_path, "_novoalign.log")
	count_df = compile_data(count_path, "_read_count.tsv")
	combined_df = pd.concat([align_df, count_df], axis=1, sort=True, join="inner")
	combined_df.to_csv(output_csv, columns=ALIGN_VARS+COUNT_VARS[1:], index_label="Sample")


if __name__ == "__main__":
	main(sys.argv)

