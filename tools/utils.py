#!/usr/bin/python
import numpy as np
from itertools import combinations
import yaml


def decompose_status2bit(n):
	"""
	Decompose the bit status 
	"""
	if n == 0:
		return None
	decomp = []
	for i in np.arange(np.floor(np.log2(n)),-1,-1):
		if n - 2**i >= 0:
			decomp.append(i)
			n -= 2**i
	return decomp


def make_combinations(lst):
	"""
	Make all possible replicate combinations
	"""
	if len(lst) < 2:
		return [lst]
	combo = []
	for i in range(len(lst), 1, -1):
		for s in combinations(lst, i):
			combo.append(s)
	return combo


def load_config(json_file):
	"""
	Load configuration file (JSON) for QC thresholding and scoring
	"""
	with open(json_file) as json_data:
		d = yaml.safe_load(json_data)
	return d


