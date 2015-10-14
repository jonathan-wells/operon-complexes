#!/usr/bin/env python3

with open('../figures/data/dataset3.csv') as infile1:
	data1 = [line.strip() for line in infile1.readlines()]	

with open('../figures/data/dataset2.txt') as infile2:
	data2 = [line.strip().split('\t') for line in infile2.readlines()]
	# for line in data2:
	# 	print(len(line))

oplen = {tuple(sorted(line[2:4])):line[-1] for line in data2[1:] if len(line) == 11}
new_data = [data1[0] + ',operon.length']
for line in data1[1:]:
	pair = tuple(sorted(line.split(',')[3:5]))
	if pair in oplen:
		line = line + ',' + oplen[pair]
	else:
		line += ',NA'
	line = line.replace('NA', '')
	new_data.append(line)
data = '\n'.join(new_data)
print(data)