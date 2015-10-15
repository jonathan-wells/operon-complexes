#!/usr/bin/env python3

import sys
sys.path.append('..')
import operon_assembly as opa

old_file = 'data/pairs_w_dimers_corrected.txt'
correct_file = 'data/dataset1.csv'

with open(correct_file) as infile:
	data = [line.strip().split(',') for line in infile.readlines()]
for line in data:
	print(line[-1])


# test = opa.intervening_genes.RandomOperons(old_file, oplevel=True)
# obs = test.calc_score()
# print(obs)
# wincount = 0
# n = 1000
# for i in range(n):
#     test.shuffle_operons()
#     ex = test.calc_score()
#     if ex <= obs:
#         wincount += 1
# print(wincount)
# print(wincount/n)