#!/usr/bin/env python3

import abundance.coabundance as co
import location.subcell as loc
import matplotlib.pyplot as plt
import random

def load_utr_sequences(UTR_type):
    """Returns dictionary of ENSG ids and corresponding UTR sequences."""
    seq_dict = {}
    filename = 'data/mrna_regulation/ensg_human_' + str(UTR_type) + 'UTR.fasta'
    with open(filename) as seqfile:
        for line in seqfile:
            if '>' in line:
                gene = line.split(':')[-1].strip()
            elif gene not in seq_dict:
                seq_dict[gene] = [line.strip()]
            else:
                seq_dict[gene].append(line.strip())
    for gene in seq_dict:
        # No good, uses longest UTR as proxy for all - prob. not valid 
        seq_dict[gene] = max(seq_dict[gene])
    return seq_dict

class GeneIDs(object):
    def __init__(self, gene_id_file):
        with open(gene_id_file) as infile:
            data = infile.readlines().split('\t')
        self._hgnc = {line[0]: line[1:] for line in data}
        self._ensp = {}
        self._ensg = {}
        self._uniprot = {}

class FeatureTable(object):
    def __init__(self, pair_file, pair_file_2, pair_file_3):
        self.header = None
        self.table = self._load_table(pair_file)
        self.strucs = [line[0] for line in self.table]
        self._add_secondary_ids(pair_file_2, pair_file_3)
        self._list_ids()

    def _load_table(self, pair_file):
        """Rough process to remove pairs from structures with >2 diff subs."""
        self.header = ['struc', 'ensp1', 'ensp2', 'type', 'interface']
        strucs = []
        data = []
        with open(pair_file) as infile:
            for line in infile:
                sline = line.split()
                if len(sline) == 5:
                    line = [sline[2], sline[0], sline[1], sline[4], sline[3]]
                    if float(sline[3]) > 0.1:
                        data.append(line)
        return data

    def filter_2subs(self):
        with open('data/pair_datasets/dimers.txt') as infile:
            dimers = [line.lower().strip() for line in infile]
        table = []
        for line in self.table:
            if line[0].split('_')[0] in dimers:
                table.append(line)
        self.table = table

    def _add_secondary_ids(self, pair_file_2, pair_file_3):
        id_head = 'ensg1 ensg2 upr1 upr2'.split()
        self.header = self.header[:3] + id_head + self.header[3:]
        upr_dict = {}
        ensg_dict = {}
        with open(pair_file_2) as infile:
            for line in infile:
                line = line.split()
                struc = line[2]
                if struc in self.strucs:
                    upr_dict[struc] = line[:2]
        with open(pair_file_3) as infile:
            for line in infile:
                line = line.split()
                struc = line[2]
                if struc in self.strucs:
                    ensg_dict[struc] = line[:2]
        table = []
        for line in self.table:
            if line[0] in upr_dict:
                line = line[:3] + upr_dict[line[0]] + line[3:]
            else:
                line = line[:3] + ['NA', 'NA'] + line[3:]
            if line[0] in ensg_dict:
                line = line[:3] + ensg_dict[line[0]] + line[3:]
            else:
                line = line[:3] + ['NA', 'NA'] + line[3:]
            table.append(line)
        self.table = table

    def _list_ids(self):
        ensps = []
        uprids = []
        ensgs = []
        for line in self.table:
            ensps += [line[1], line[2]]
            uprids += [line[5], line[6]]
            ensgs += [line[3], line[4]]
        self.ensps = sorted(set(ensps))
        self.ensgs = sorted(set(ensgs))
        self.uprids = sorted(set(uprids))

    def add_expression_correlations(self, ptype):
        self.header.append(ptype)
        comat = co.CorrelationMatrix(self.uprids, ptype)
        for line in self.table:
            try:
                cor = comat.calc_correlation(line[5], line[6])[0]
                line.append(str(round(cor, 2)))
            except:
                line.append('NA')

    def add_localisation_data(self):
        self.header.append('loc')
        cell1 = loc.Cell('experiment')
        cell2 = loc.Cell('knowledge')
        for line in self.table:
            try:
                ploc = cell1.pairwise_location(line[1], line[2], 3)
                line.append(str(round(ploc, 2)))
            except:
                try:
                    ploc = cell2.pairwise_location(line[1], line[2], 3)
                    line.append(str(round(ploc, 2)))
                except:
                    line.append('NA')

    def add_utr_data(self, utr_type):
        self.header.append(str(utr_type)+'utr')
        utrs = load_utr_sequences(utr_type)
        for line in self.table:
            try:
                ensp1 = len(utrs[line[3]])
                ensp2 = len(utrs[line[4]])
                avg = (ensp1 + ensp2)/2
                diff = ((ensp1 - ensp2)**2)**0.5
                line.append(str(diff))
            except:
                line.append('NA')

    def add_tf_reg_data(self):
        self.header.append('tfreg')
        pass

    def randomise_pairs(self):
        self.header = 'struc ensp1 ensp2 ensg1 ensg2 upr1 upr2'.split()
        p1_list = [(i[1], i[3], i[5]) for i in self.table]
        p2_list = [(i[2], i[4], i[6]) for i in self.table]
        table = []
        for line in self.table:
            p1 = random.choice(p1_list)
            p2 = random.choice(p2_list)
            table.append(['rand', p1[0], p2[0], p1[1], p2[1], p1[2], p2[2]])
        self.table = table

    def __str__(self):
        header = '\t'.join(self.header)
        ptable = ['\t'.join(line) for line in self.table]
        ptable = '\n'.join(sorted(ptable))
        return header + '\n' + ptable


def main():
    ft = FeatureTable('data/pair_datasets/ensp_human_pairs.txt',
                      'data/pair_datasets/upr_human_pairs.txt',
                      'data/pair_datasets/ensg_human_pairs.txt')
    # ft.filter_2subs()
    # ft.randomise_pairs()
    ft.add_expression_correlations('mRNA')
    ft.add_expression_correlations('protein')
    ft.add_localisation_data()
    ft.add_utr_data(3)
    ft.add_utr_data(5)
    print(ft)

if __name__ == '__main__':
   main()