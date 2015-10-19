#!/usr/bin/env python3

import sys
sys.path.append('..')
from operon_assembly import random_operons as ro

class SignificanceTests(ro.RandomOperons):
	"""Provides methods to test the significance of 
	"""

	_strucs = 'data/dataset1.csv'
	_yeast2hybrid = 'data/ecoli_y2h_pairwise.txt'
	_struc_indices = [3, 4, 6, 14, 15, 16, 18]
	_y2h_indices = [0, 1, 6, 4, 2, 3, 5]
	_struc_sep = ','
	_y2h_sep = '\t'
	
	def __init__(self, dtype):
		"""Initiates data based on either y2h or structural data."""
		if dtype == 'struc':
			data = self._strucs
			inds = self._struc_indices
			sep = self._struc_sep
		elif dtype == 'y2h':
			data = self._yeast2hybrid
			inds = self._y2h_indices
			sep = self._y2h_sep
		else:
			raise ValueError('use "struc" or "y2h"')
		super().__init__(data, inds, sep, oplevel=True)
		self.dtype = dtype

	def fractional_test(self, trials):
		"""Tests significance of observing fraction of interacting genes that
		are also adjacent.

		Returns:
			observed fraction of interacting genes that are adjacent
			mean expected fraction
			p-value
		"""
		observed = self.test.calc_fraction()
		successes = 0
		mean = 0
		for i in range(trials):
			self.test.shuffle_operons()
			expected = self.test.calc_fraction()
			mean += expected
			if expected >= observed:
				successes += 1
		observed = round(observed, 2)
		mean_expected = round(mean/trials, 2)
		pval = round(successes/trials, 2)
		self.__init__(self.dtype)
		return (observed, mean_expected, pval)

	def total_intervening_test(self, trials):
		"""Tests significance of observing at most x number of intervening
		genes between all interacting pairs in dataset.

		Returns:
			observed number of intervening genes between interacting pairs.
			mean expected number based on random positioning within operon.
			p-value
		"""
		observed = self.calc_intervening()
		successes = 0
		mean = 0
		for i in range(trials):
			self.shuffle_operons()
			expected = self.calc_intervening()
			mean += expected
			if expected <= observed:
				successes += 1
		observed = round(observed, 2)
		mean_expected = round(mean/trials, 2)
		pval = round(successes/trials, 2)
		self.__init__(self.dtype)
		return (observed, mean_expected, pval)

	def intervening_distribution(self, trials):
		"""See total_intervening_test(). Returns raw values of that."""
		observed = self.calc_intervening(printr=True)[1]
		expected = []
		for i in range(trials):
			self.shuffle_operons()
			expected += self.calc_intervening(printr=True)[1]
		return observed, expected

	def summary(self):
		"""Summarises information about operons in data"""
		num_ops = len(self.oplen.keys())
		total_oplen = sum(self.oplen.values())
		print(num_ops, total_oplen, total_oplen/num_ops)


if __name__ == '__main__':
	struc = SignificanceTests('struc')
	# complex 1 operons: '116617', '580994'
	struc.filter_operons('116617', '580994', exclude=False)
	struc.summary()
	# print(struc)