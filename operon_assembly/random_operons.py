#!/usr/bin/env python3

import random
import itertools
import math

"""Think I used this to calculate some sort of null model with shuffled ops?"""

class RandomOperons(object):
    def __init__(self, filename, data_indices, sep,
        complevel=False, oplevel=False):
        """Loads information about operons.
        Arguments:
            filename - input file containing operon info
            data_indices - [g1, g2, int, opID, g1pos, g2pos, oplength]
                where gxpos is the position of gene x in operon
            sep - separator to use when splitting lines.
            complevel - shuffle at level of complexes within operons.
            oplevel - shuffle at level of complete operons.
        """
        self.levels = complevel, oplevel
        self.ind = data_indices
        with open(filename) as infile:
            self.data = [l.strip().split(sep) for l in infile.readlines()[1:]]
            self.data = [l for l in self.data if l[-1] != '']
        self._init_int_dict()
        self._init_op_genes()
        self._init_op_len()

    def _init_op_len(self):
        """Dictionary of operons and their lengths."""
        self.oplen = {}
        for line in self.data:
            self.oplen[line[self.ind[3]]] = int(line[self.ind[6]])

    def _init_int_dict(self):
        """Dictionary of operons. Each operon points to a dictionary of gene
        pairs. Each gene pair returns the size of the interface between those
        genes.
        """
        self.int_dict = {}
        for line in self.data:
            op = line[self.ind[3]]
            pair = tuple(sorted([line[self.ind[0]], line[self.ind[1]]]))
            interface = line[self.ind[2]]
            # Allows use of ecoli Y2H data
            if interface == 'TRUE':
                interface = 1000.0
            elif interface == 'FALSE':
                interface = 0.0
            else:
                interface = float(interface)
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
            op = line[self.ind[3]]
            g1 = line[self.ind[0]]
            g2 = line[self.ind[1]]
            pos1 = int(line[self.ind[4]])
            pos2 = int(line[self.ind[5]])
            if op not in self.op_genes:
                self.op_genes[op] = {g1: pos1, g2: pos2}
            else:
                self.op_genes[op][g1] = pos1
                self.op_genes[op][g2] = pos2

    def filter_operons(self, *args, exclude=True):
        """Bit messy but okay"""
        for op in args:
            if op not in self.op_genes:
                continue
            if exclude == True:
                self.op_genes.pop(op)
                self.int_dict.pop(op)
                self.oplen.pop(op)
            # Limits operons exclusively to the ones supplied in args
            elif exclude == False:
                for op in list(self.op_genes):
                    if op not in args:
                        self.op_genes.pop(op)
                        self.int_dict.pop(op)
                        self.oplen.pop(op)


    def calc_intervening(self, printr=False):
        """Calculate total number of intervening genes separating interacting
        pairs across dataset.
        """
        score = 0
        intervenings = []
        for op in self.int_dict:
            for pair in self.int_dict[op]:
                genes = tuple(pair)
                if self.int_dict[op][pair] > 200:  # i.e. pair form interface
                    pos1 = self.op_genes[op][genes[0]]
                    pos2 = self.op_genes[op][genes[1]]
                    intervening = math.sqrt((pos1 - pos2)**2) - 1
                    score += intervening
                    intervenings.append(intervening)
        if printr == True:
            return score, intervenings
        else:
            return score

    def calc_fraction(self):
        """Calculate total fraction of interacting pairs that are adjacent."""
        adj_score = 0
        int_score = 0
        for op in self.int_dict:
            for pair in self.int_dict[op]:
                genes = tuple(pair)
                if self.int_dict[op][pair] > 200:  # i.e. genes form interface
                    pos1 = self.op_genes[op][genes[0]]
                    pos2 = self.op_genes[op][genes[1]]
                    intervening = math.sqrt((pos1 - pos2)**2) - 1
                    if intervening == 0:  # i.e. genes adjacent within operon
                        adj_score += 1
                    int_score += 1
        return adj_score/int_score

    def shuffle_operons(self):
        complevel = self.levels[0]
        oplevel = self.levels[1]
        for op in self.op_genes:
            if complevel == True:
                pos_list = list(self.op_genes[op].values())
            elif oplevel == True:
                pos_list = [i for i in range(1, (self.oplen[op] + 1))]
            random.shuffle(pos_list)
            index = 0
            for gene in self.op_genes[op]:
                self.op_genes[op][gene] = pos_list[index]
                index += 1

    def find_multicomp_operons(self):
        """Prints out any operons that contain multiple PDB ids. Could be due
        to either multiple real complexes, or just redundant PDB ids.
        """
        self.ocomp = {}
        for line in self.data:
            if line[14] not in self.ocomp:
                self.ocomp[line[14]] = [line[0]]
            else:
                self.ocomp[line[14]].append(line[0])
        for comp in self.ocomp:
            self.ocomp[comp] = set(self.ocomp[comp])
            if len(self.ocomp[comp]) != 1:
                print(comp, self.ocomp[comp])

    def __repr__(self):
        return(str('\n'.join(['\t'.join(line) for line in self.data])))


if __name__ == '__main__':
    test = RandomOperons()
