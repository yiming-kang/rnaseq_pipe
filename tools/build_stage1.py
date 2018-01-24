#!/usr/bin/python
import sys
import argparse
import pandas as pd


def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--samples', 
						help='Sample summary metadata file.')
	parser.add_argument('-i', '--genome_index', 
						help='Genome index file built for the aligner, e.g. *.nix for novoalign.')
	parser.add_argument('-r', '--reference_gtf', 
						help='Annotation reference file in GTF/GFF3 format.')
	parser.add_argument('-g', '--gene_list',
						help='Gene list.')
	return parser.parse_args(argv[1:])


def build_alignment(samples, genome_index):
	align_dict = {}
	for i,row in samples.iterrows():
		sample_id = row['GENOTYPE'] +'-'+ row['SAMPLE']
		prereq = row['FILE']
		target_dir = 'alignment/novoalign/'+ sample_id
		target_file = target_dir +'/aligned_reads_sorted.bam'
		recipe = 'rm -rf '+ target_dir +'; mkdir '+ target_dir +'; '
		recipe += 'novoalign -o SAM -d '+ genome_index +' -f '+ prereq +' 2> '+ target_dir +'/novoalign.log | samtools view -bS > '+ target_dir +'/aligned_reads.bam; '
		recipe += 'novosort '+ target_dir +'/aligned_reads.bam > '+ target_dir +'/aligned_reads_sorted.bam 2> '+ target_dir +'/novosort.log; '
		recipe += 'rm '+ target_dir +'/aligned_reads.bam'
		align_dict[target_file] = {'prereq': prereq, 'recipe': recipe}
	return align_dict


def build_expression_quantification(samples, reference_gtf, gene_list):
	expr_dict = {}
	for i,row in samples.iterrows():
		sample_id = row['GENOTYPE'] +'-'+ row['SAMPLE']
		prereq = 'alignment/novoalign/'+ sample_id +'/aligned_reads_sorted.bam'
		target_dir = 'expression/stringtie/'+ sample_id
		target_file = 'expression/stringtie_fpkm/'+ sample_id +'.expr'
		recipe = 'rm -rf '+ target_dir +'; mkdir '+ target_dir +'; '
		recipe += 'stringtie '+ prereq +' -G '+ reference_gtf +' -e -o '+ target_dir +'/stringtie_out.gtf -A '+ target_dir +'/gene_abundances.tab; '
		recipe += 'python tools/stringtie2fpkm.py '+ target_dir +'/gene_abundances.tab '+ gene_list +' > '+ target_file
		expr_dict[target_file] = {'prereq': prereq, 'recipe': recipe}
	return expr_dict


def main(argv):
	parsed = parse_args(argv)
	samples = pd.read_csv(parsed.samples, delimiter='\t', dtype=str)
	
	align_dict = build_alignment(samples, parsed.genome_index)
	expr_dict = build_expression_quantification(samples, parsed.reference_gtf, parsed.gene_list)
	
	for target in sorted(align_dict.keys()):
		sys.stdout.write('%s: %s\n\t%s\n' % (target, align_dict[target]['prereq'], align_dict[target]['recipe']))
	sys.stdout.write('NOVOALIGN_ALIGNMENTS = %s\n' % ' '.join(sorted(align_dict.keys())))
	for target in sorted(expr_dict.keys()):
		sys.stdout.write('%s: %s\n\t%s\n' % (target, expr_dict[target]['prereq'], expr_dict[target]['recipe']))
	sys.stdout.write('STRINGTIE_EXPRESSIONS = %s\n' % ' '.join(sorted(expr_dict.keys())))
	sys.stdout.write('all: $(NOVOALIGN_ALIGNMENTS) $(STRINGTIE_EXPRESSIONS)')

if __name__ == '__main__':
	main(sys.argv)
