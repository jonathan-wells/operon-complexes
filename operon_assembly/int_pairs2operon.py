#!/usr/bin/env python3

import re
import operator
import ast
import collections
import glob
from itertools import combinations
import math

class IntPairs2Operon(object):

    """Methods to produce table of interacting gene pairs and operon info."""

    def __init__(self, opr_filename, gene_pairs, out_file, tax_id):
        self.opr_filename = opr_filename
        self.gene_pairs = gene_pairs
        self.out_file = out_file
        self.tax_id = tax_id
        self._compile_op_dict()

    def _compile_op_dict(self):
        """Generate dictionaries of operons containg gene order information."""
        with open(self.opr_filename) as file:
            data = file.readlines()
        # Pattern matches all the info about each gene in .opr file
        gene_pattern = re.compile(r'(\d*)\s*(\S*)\s*(\S*)\s*(\d*)\s*(\d*)'
                                  '\s*([+-])\s*(\d*)\s*(\S *)\s*(.*)')
        self.op_dict = {}
        error_list = []
        for line in data:  # Formerly [1:]
            if line == '\n':
                continue
            else:
                try:
                    hit = gene_pattern.match(line)
                    operon_id = hit.group(1)
                    gene_sym = hit.group(2)
                    strand = hit.group(6)
                except:
                    error_list.append(line)
                    continue
            if operon_id not in self.op_dict:
                self.op_dict[operon_id] = [[gene_sym], strand] # Note strand
            else:
                self.op_dict[operon_id][0].append(gene_sym)
        # Strips single gene operons and reverses anti-sense gene orders.
        self.op_dict_final = collections.OrderedDict()
        for ident in self.op_dict:
            if self.op_dict[ident][1] == '-':
                self.op_dict[ident][0].reverse()
            if len(self.op_dict[ident][0]) > 1:
                self.op_dict_final[ident] = self.op_dict[ident][0]
        self.error_list = error_list
        self.op_dict = self.op_dict_final
        # print(error_list)
        # return self.op_dict_final

    def pos_finder(self, genA, genB, gene_list):
        """Find the coordinates of geneA and B in the operon."""
        posA = gene_list.index(genA)+1
        posB = gene_list.index(genB)+1
        coords = [(posA, posB)]
        if posA < posB:
            coords.append(genA)
        else:
            coords.append(genB)
        return coords

    def check_gene_pairs(self):
        """For a given operon dictionary checks to see whether interacting
        gene pairs are present.
        """
        with open(self.gene_pairs) as file:
            gene_data = file.readlines()
        # lin_pattern = re.compile(r'(\w*)\s([\w]*)\s([A-Z0-9]*)\s([A-Z0-9]*)')
        structure_dict = {}
        results = []
        for line in gene_data:
            line = line.split()
            try:
                # hit = lin_pattern.match(line)
                struc = line[0].split('_')[0]
                org = 'NA'
                genA = line[2]
                genB = line[3]
                # print(genA, genB)
            except:
                continue
            # Beauty
            for op in self.op_dict:
                if genA in self.op_dict[op] and genB in self.op_dict[op]:
                    coords = self.pos_finder(genA, genB, self.op_dict[op])
                    self.out_file.write("{:<8}".format(struc))
                    self.out_file.write("{:<11}".format(op))
                    self.out_file.write("{:<10}".format(genA))
                    self.out_file.write("{:<10}".format(genB))
                    self.out_file.write("{:<10}".format(str(coords[0])))
                    self.out_file.write("{:<10}".format(str(coords[1])))
                    self.out_file.write("{:<8}".format(self.tax_id))
                    self.out_file.write(org+'\n')
                    

def generate_results_file():
    gene_pairs = 'data/prokaryotic_gene_pairs/pairs_missing.txt'
    with open('data/test.txt') as file:
        data = file.readlines()
    x = ast.literal_eval(data[0])
    with open('data/int_pair2operon/test.txt', 'a') as out_file:
        header_positions = "{:<8}{:<11}{:<10}{:<10}{:<10}{:<10}{:<8}"
        out_file.write(header_positions.format('Struc', 'Opr_ID', 'Gene_A',
                                               'Gene_B', 'Coords', 'Tr. 1st',
                                               'Tax_ID'))
        out_file.write('Organism\n')
        bac = glob.glob('./data/operon_datasets/bac_uniprot_opr_files/*')
        arc = glob.glob('./data/operon_datasets/arc_uniprot_opr_files/*')
        ext = glob.glob('./data/operon_datasets/ext_uniprot_opr_files/*')
        mis = glob.glob('./data/operon_datasets/missing_opr/*')
        bac_f = [f.split('/')[-1] for f in bac]
        arc_f = [f.split('/')[-1] for f in arc]
        ext_f = [f.split('/')[-1] for f in ext]
        mis_f = [f.split('/')[-1] for f in mis]
        for item in x:
            f = 'uniprot_'+item[0]+'a.opr'
            if f in bac_f:
                location = bac[bac_f.index(f)].strip('./')
            elif f in arc_f:
                location = arc[arc_f.index(f)].strip('./')
            elif f in ext_f:
                location = ext[ext_f.index(f)].strip('./')
            elif f in mis_f:
                location = mis[mis_f.index(f)].strip('./')
            final_results = IntPairs2Operon(location, gene_pairs,
                                            out_file, item[2])
            final_results.check_gene_pairs()
        for item in x:
            f = 'uniprot_'+item[0]+'b.opr'
            if f in bac_f:
                location = bac[bac_f.index(f)].strip('./')
            elif f in arc_f:
                location = arc[arc_f.index(f)].strip('./')
            elif f in ext_f:
                location = ext[ext_f.index(f)].strip('./')
            elif f in mis_f:
                location = mis[mis_f.index(f)].strip('./')
            final_results = IntPairs2Operon(location, gene_pairs,
                                            out_file, item[2])
            x = final_results.check_gene_pairs()
    return

def sort_results_file(results_file):
    """Sorts file by tax_id so that results include more detailed strain info"""
    file = open(results_file)
    data = file.readlines()
    sorted_copy = []
    for line in data[1:]:
        sorted_copy.append(line.split())
    for line in sorted_copy[1:]:
        line[7] = int(line[7])
    # UGLY, i think this can be done with less code - do you have to reprint each .split()?
    with open('data/int_pair2operon/missing_sorted.txt','w') as out_file:
        out_file.write(data[0])
        for line in sorted(sorted_copy[1:], key=operator.itemgetter(7, 0), 
                           reverse=True):
            out_file.write("{:<8}".format(line[0]))
            out_file.write("{:<11}".format(line[1]))
            out_file.write("{:<10}".format(line[2]))
            out_file.write("{:<10}".format(line[3]))
            out_file.write("{:<10}".format(line[4]+' '+line[5]))
            out_file.write("{:<10}".format(line[6]))
            out_file.write("{:<8}".format(line[7]))
            out_file.write(line[8]+'\n')
    file.close()
    return
    
def clean_results_file(results_file):
    """Removes duplicate entries where gene_pairs have been found in multiple 
    .opr files."""
    file = open(results_file)
    data = file.readlines()
    clean_copy = []
    gpt_list = []
    with open('data/int_pair2operon/missing_clean.txt','w') as outfile:
        for line in data:
            contents = line.split()
            gene_pair = (contents[2], contents[3])
            if gene_pair not in gpt_list:
                gpt_list.append(gene_pair)
                clean_copy.append(line)
        for item in clean_copy:
           outfile.write(item) # +''
    file.close()
    
def reformat_in_exp_order(filename):
    """Reformats results in order of gene expression."""
    with open(filename) as file:
        data = file.readlines()
    with open('data/int_pair2operon/missing_final.txt','w') as outfile:
        header = '{:<8}{:<10}{:<10}{:<10}{:<10}{:<10}{}'
        outfile.write(header.format('Struc', 'Opr_ID', 'Gene_A', 'Gene_B', 
                                    'Coords', 'Tax_ID', 'Organism')+'\n')
        for line in data[1:]:
            sline = line.split('  ')
            sline = list(filter(None, sline))   
            sline = [x.strip() for x in sline]    
            # i.e if geneA expressed after geneA, flip stuff about
            pos = ast.literal_eval(sline[4])
            genA = sline[2]; genB = sline[3]
            if pos[0] > pos[1]:
                sline[2] = genB; sline[3] = genA
                sline[4] = (pos[1], pos[0])     
            s = sline
            # Pretty prints everything
            if len(s) == 8:
                form = '{:<8}{:<10}{:<10}{:<10}{:<10}{:<10}{}'
                line = form.format(s[0], s[1], s[2], s[3], str(s[4]), s[6], s[7])
            outfile.write(line+'\n')


class EcoliY2H(object):
    """For looking at binarpy PPIs, non-structural."""
    def __init__(self, y2h_file, op_file_id):
        self._init_ppi_pairs(y2h_file)
        self._init_op_dict(op_file_id)
        
    def _init_ppi_pairs(self, y2h_file):
        with open(y2h_file) as file:
            data = file.readlines()
        pair_list = []
        for line in data[1:]:
            line = line.split()
            if 'uniprot' in line[0] and 'uniprot' in line[1]:
                geneA = line[0].split('uniprotkb:')[1]
                geneB = line[1].split('uniprotkb:')[1]
                if geneA != geneB:
                    pair_list.append(tuple(sorted([geneA, geneB])))
        self.ppi_pairs = set(pair_list)
        
    def _init_op_dict(self, fid):
        uniprot_files = ['uniprot_'+fid+'a.opr', 'uniprot_'+fid+'b.opr']
        final_op_dict = {}
        for f in uniprot_files:
            fname = 'data/operon_datasets/bac_uniprot_opr_files/'+f
            I2O = IntPairs2Operon(fname, None, None, None)
            for op in I2O.op_dict:
                if op not in final_op_dict:
                    final_op_dict[op] = I2O.op_dict[op]
        for op in list(final_op_dict):
            # if len(final_op_dict[op]) != len(set(final_op_dict[op])):
            #     final_op_dict.pop(op)
            #     continue
            for prot in final_op_dict[op]:
                if len(prot) != 6:
                    final_op_dict.pop(op)
                    break
        self.op_dict = final_op_dict
    
    def get_oplen(self, op):
        return len(self.op_dict[op])
    
    def pos_finder(self, genA, genB, gene_list):
        """Find the coordinates of geneA and B in the operon."""
        posA = gene_list.index(genA)+1
        posB = gene_list.index(genB)+1
        coords = [(posA, posB)]
        if posA < posB:
            coords.append(genA)
        else:
            coords.append(genB)
        return coords
    
    def get_potential_operon_ppis(self):
        print('geneA_id', 'geneB_id', 'geneA_position', 'geneB_position',
              'intervening_genes', 'operon_id', 'operon_length',
              'ppi_detected', sep='\t')
        for op in self.op_dict:
            oplen = self.get_oplen(op)
            op_combinations = list(combinations(self.op_dict[op], 2))
            op_combinations = [tuple(sorted(i)) for i in op_combinations]
            for pair in op_combinations:
                if pair[0] == pair[1]:
                    continue
                coords = self.pos_finder(pair[0], pair[1], self.op_dict[op])
                inter = str(int(math.sqrt((coords[0][0] - coords[0][1])**2) - 1))
                coords = tuple(str(i) for i in coords[0])
                line = ['\t'.join(pair), '\t'.join(coords), inter, op, str(oplen)]
                if pair in self.ppi_pairs:
                    line.append('TRUE')
                else:
                    line.append('FALSE')
                print('\t'.join(line))
            
            
    
if __name__ == '__main__':
    eY2H = EcoliY2H('data/ecoli_y2h_raw.txt', '1034')
    eY2H.get_potential_operon_ppis()
    