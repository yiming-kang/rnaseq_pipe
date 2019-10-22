import sys
import os
import argparse
from glob import glob


FASTQ_TYPES = ["fastq.gz", "fq.gz"]
FEATURE_TYPE_DICT = {"gff": "gene", "gtf": "CDS"}


def parse_args(argv):
	parser = argparse.ArgumentParser(description="This script generates sbatch script and submits sbatch job.")
	parser.add_argument("-f", "--fastq_path", required=True,
						help="[Required] Directory path of fastq files of type: {}".format(FASTQ_TYPES))
	parser.add_argument("-i", "--genome_index_file", required=True,
						help="[Required] File path of genome index. This is specific to the aligner.")
	parser.add_argument("-s", "--strandness", required=True,
						help="[Required] Specify 'yes', 'no', or 'reverse'. For NEB kit, use 'reverse'.")
	parser.add_argument("--gene_annotation_file", default=None,
						help="[Optional]  File path of gene annotation. By default (if not specified), it will look for .gff or .gtf file in the same directory and has same filename as genome index file.")
	parser.add_argument("--annotation_feature_type", default=None,
						help="[Optional]  Feature type to use for reads counting. By default (if not specified), it will use this dictionary {} based on the annotation file type.".format(FEATURE_TYPE_DICT))
	parser.add_argument("--output_path", default=None,
						help="[Optional] Directory path of output data. By default (if not specified), it will write output files in the same directory as the input fastq files.")
	parser.add_argument("--user_email", default=None,
						help="[Optional] Email for job status notification.")
	parser.add_argument("--align_only", action="store_true",
						help="[Optional] Set this flag to only align reads.")
	args = parser.parse_args(argv[1:])
	return args


def find_annotation_file(path_prefix):
	parsed = None
	for suffix in FEATURE_TYPE_DICT.keys():
		file_path = path_prefix + "." + suffix
		if os.path.exists(file_path):
			parsed = file_path
	if parsed is None:
		sys.exit("Error: annotation file not found.")
	return parsed


def get_feature_type(file_path):
	suffix = file_path.split(".")[-1]
	if suffix not in FEATURE_TYPE_DICT.keys():
		sys.exit("Error: annotation file has incorrect suffix.")
	return FEATURE_TYPE_DICT[suffix]


def write_fastq_list(dir_path, out_file):
	file_paths = []
	for suffix in FASTQ_TYPES:
		file_paths += glob(dir_path + "/*." + suffix)
	with open(out_file, "w") as f:
		for file_path in file_paths:
			f.write("{}\n".format(file_path))
	return len(file_paths)


def write_job_script(job_file, output_path, fastq_list_file, num_fastqs, geno_idx_file, gene_ann_file, feat_type, strandness, align_only):
	with open(job_file, "w") as f:
		f.write("#!/bin/bash\n")
		f.write("#SBATCH -N 1\n")
		f.write("#SBATCH --cpus-per-task=8\n")
		f.write("#SBATCH --mem=12G\n")
		f.write("#SBATCH --array=1-{0}%{1}\n".format(num_fastqs, min(num_fastqs, 50)))
		f.write("#SBATCH -D ./\n")
		f.write("#SBATCH -o log/mblab_rnaseq_%A_%a.out\n")
		f.write("#SBATCH -e log/mblab_rnaseq_%A_%a.err\n")
		f.write("#SBATCH -J mblab_rnaseq\n\n")
		f.write("ml novoalign/3.07.00\n")
		f.write("ml samtools/1.6\n")
		f.write("ml htseq/0.9.1\n")
		f.write("read fastq_file < <( sed -n ${{SLURM_ARRAY_TASK_ID}}p {} ); set -e\n\n".format(fastq_list_file))
		f.write("mkdir -p {}\n".format(output_path))
		f.write("sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; novoalign -c 8 -o SAM -d {0} -f ${{fastq_file}} 2> log/${{sample}}_novoalign.log | samtools view -bS > {1}/${{sample}}_aligned_reads.bam\n".format(geno_idx_file, output_path))
		f.write("sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; novosort --threads 8 {0}/${{sample}}_aligned_reads.bam > {0}/${{sample}}_sorted_aligned_reads.bam 2> log/${{sample}}_novosort.log\n".format(output_path))
		if align_only is False:
			f.write("sample=${{fastq_file##*/}}; sample=${{sample%.f*q.gz}}; htseq-count -f bam -i ID -s {1} -t {2} {0}/${{sample}}_sorted_aligned_reads.bam {3} > {0}/${{sample}}_read_count.tsv 2> log/${{sample}}_htseq.log\n".format(output_path, strandness, feat_type, gene_ann_file))


def main(argv):
	args = parse_args(argv)
	fastq_path = args.fastq_path
	output_path = args.output_path 
	geno_idx_file = args.genome_index_file
	gene_ann_file = args.gene_annotation_file
	ann_feat_type = args.annotation_feature_type
	strandness = args.strandness
	user_email = args.user_email
	align_only = args.align_only

	## Parse default variables
	if output_path is None:
		output_path = fastq_path

	if gene_ann_file is None:
		gene_ann_prefix = ".".join(geno_idx_file.split(".")[:-1])
		gene_ann_file = find_annotation_file(gene_ann_prefix)
	 
	if ann_feat_type is None:
		feat_type = get_feature_type(gene_ann_file)

	## Write sbatch script
	fastq_list_file = "fastq_list.txt"
	sbatch_job_file = "mblab_rnaseq.sbatch"
	os.system("mkdir -p log/")
	num_fastqs = write_fastq_list(fastq_path, fastq_list_file)
	write_job_script(sbatch_job_file, output_path,
					fastq_list_file, num_fastqs, 
					geno_idx_file, gene_ann_file, feat_type, 
					strandness, align_only)

	## Submit sbatch job
	if user_email is None:
		os.system("sbatch {}".format(sbatch_job_file))
	else:
		os.system("sbatch --mail-type=END,FAIL --mail-user={0} {1}".format(user_email, sbatch_job_file))

if __name__ == "__main__":
	main(sys.argv)

