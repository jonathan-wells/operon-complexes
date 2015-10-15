#!/usr/bin/env python3

import random
import itertools
import math

"""Think I used this to calculate some sort of null model with shuffled ops?"""

class RandomOperons(object):
    def __init__(self, filename, complevel=False, oplevel=False):
        self.levels = complevel, oplevel
        with open(filename) as file:
            self.data = file.readlines()[1:]
        self._init_int_dict()
        self._init_op_genes()
        self._init_op_len()
        self._filter_operons('116617')
        
    def _init_op_len(self):
        """Dictionary of operons and their lengths."""
        self.oplen = {}
        for line in self.data:
            line = line.split()
            if len(line) > 8:
                self.oplen[line[8]] = int(line[11])
            
    def _init_int_dict(self):
        """Dictionary of operons. Each operon points to a dictionary of gene 
        pairs. Each gene pair returns the size of the interface between those 
        genes.
        """
        self.int_dict = {}
        for line in self.data:
            line = line.split()
            if len(line) > 8:
                op = line[8]
                pair = tuple(sorted([line[2], line[3]]))
                interface = float(line[5])
                if op not in self.int_dict:
                    self.int_dict[op] = {pair: interface}
                else:
                    self.int_dict[op][pair] = interface
    
    def _init_op_genes(self):
        """Dictionary of operons. Each operon points to a dictionary of genes.
        Each gene returns its position in the operon.
        """
        self.op_genes = {}
        for line in self.data:
            line = line.split()
            if len(line) > 8:
                op = line[8]
                g1 = line[2]
                g2 = line[3]
                pos1 = int(line[9].strip('(,'))
                pos2 = int(line[10].strip(')'))
                if op not in self.op_genes:
                    self.op_genes[op] = {g1: pos1, g2: pos2}
                else:
                    self.op_genes[op][g1] = pos1
                    self.op_genes[op][g2] = pos2
                
    def _filter_operons(self, *args):
        """Bit messy but okay"""
        for op in args:
            self.op_genes.pop(op)
            self.int_dict.pop(op)
        if self.levels[0] == True or self.levels[1] == True:
            for op in list(self.op_genes):
                adcount = 0
                nonadcount = 0
                sortedpos = sorted(self.op_genes[op].values())
                for i in range(1,len(sortedpos)):
                    if sortedpos[i] - sortedpos[i-1] == 1:
                        adcount += 1
                        break
                if sortedpos[-1] - sortedpos[0] != 1:
                    nonadcount += 1
                if adcount == 0 or nonadcount == 0:
                    self.op_genes.pop(op)
                    self.int_dict.pop(op)
        if self.levels[1] == True:
            for op in list(self.op_genes):
                if len(self.op_genes[op]) < 3:
                    self.op_genes.pop(op)
                    self.int_dict.pop(op)
                

    def calc_score(self):
        score = 0
        for op in self.int_dict:
            for pair in self.int_dict[op]:
                genes = tuple(pair)
                if self.int_dict[op][pair] >= 500:
                    pos1 = self.op_genes[op][genes[0]]
                    pos2 = self.op_genes[op][genes[1]]
                    intervening = math.sqrt((pos1 - pos2)**2) - 1
                    score += intervening
        return score
        
    def shuffle_operons(self):
        complevel = self.levels[0]
        oplevel = self.levels[1]
        for op in self.op_genes:
            if complevel == True:
                pos_list = list(self.op_genes[op].values())
            elif oplevel == True:
                pos_list = [i for i in range(1, (self.oplen[op] + 1))]
            random.shuffle(pos_list)
            # print('')
            index = 0
            for gene in self.op_genes[op]:
                # print(gene, self.op_genes[op][gene], pos_list[index])
                self.op_genes[op][gene] = pos_list[index]
                index += 1
                        
    def find_multicomp(self):
        self.ocomp = {}
        for line in self.data:
            line = line.split()
            if len(line) > 8:
                if line[8] not in self.ocomp:
                    self.ocomp[line[8]] = [line[0].split('_')[0], line[1].split('_')[0]]
                else:
                    self.ocomp[line[8]].append(line[0].split('_')[0])
                    self.ocomp[line[8]].append(line[1].split('_')[0])
        for comp in self.ocomp:
            self.ocomp[comp] = set(self.ocomp[comp])
            if len(self.ocomp[comp]) != 1:
                print(comp, self.ocomp[comp])


def main(n, complevel=False, oplevel=False):
    filename = 'data/prokaryotic_gene_pairs/pairs_w_dimers_corrected.txt'
    test = RandomOperons(filename, complevel, oplevel)
    obs = test.calc_score()
    wincount = 0
    for i in range(n):
        test.shuffle_operons()
        ex = test.calc_score()
        if ex <= obs:
            wincount += 1
    print(wincount)
    print(wincount/n)
    
if __name__ == '__main__':
    # main(10000, complevel=True)
    main(1000000, oplevel=True)
    