#!/usr/bin/env python3

"""pretty sure this file is used to find the position of all genes in an
operon. Not 100%. so make sure and rename file to something better. also
write explanation in readme if needs be. at the very least write some
docstrings!
"""


import glob

def find_missing():
    with open('data/prokaryotic_gene_pairs/final_pairs_dataset.txt') as file:
        data = file.readlines()
    with open('data/prokaryotic_gene_pairs/bacterial_gene_pairs.txt') as file:
        bac_data = file.readlines()
    with open('data/prokaryotic_gene_pairs/archaeal_gene_pairs.txt') as file:
        arc_data = file.readlines()
    with open('data/prokaryotic_gene_pairs/extended_pairs.txt') as file:
        ext_data = file.readlines()
    with open('data/prokaryotic_gene_pairs/ribosomal_pairs.txt') as file:
        rib_data = file.readlines()
    old_data = bac_data + arc_data + ext_data + rib_data
    pairs = []
    for line in data:
        line = line.split()
        pair = ' '.join(line[2:4])
        if len(line) == 8: # I.e not in operon or not analysed
            pairs.append(pair)
    old_pairs = []
    for line in old_data:
        line = line.split()
        old_pair = ' '.join(line[2:4])
        old_pairs.append(old_pair)
    results = []
    for pair in pairs:
        rev_pair = ' '.join(pair.split()[::-1])
        if pair not in old_pairs and rev_pair not in old_pairs:
            results.append(pair)
        
    return results

def remove_missing(missing):
    with open('data/prokaryotic_gene_pairs/final_pairs_dataset.txt') as file:
        data = file.readlines()
    new_data = []
    for line in data:
        pair = ' '.join(line.split()[2:4])    
        if pair in missing:
            new_data.append(line.strip())
    for line in new_data:
        print(line)
    
def check_opr_files():
    with open('data/test.txt') as file:
        missing_data = file.readlines()
    directories = glob.glob('./data/operon_datasets/*')
    files = []
    for di in directories:
        files += glob.glob(di+'/*')
    existing = []
    for path in files:
        existing.append(path.split('/')[-1].strip('.opr'))
    present = []
    for opr in missing_data:
        if opr.strip() in existing:
            present.append((opr.strip()))
    to_download = []
    for opr in missing_data:
        if opr.strip() not in present:
            print(opr.strip())
            to_download.append((opr.strip()))
    
    

def append_new_operons(f1, f2):
    with open(f1) as f1, open(f2) as f2:
        opr_info = f1.readlines()
        complete = f2.readlines()
    opr = {}
    for i in opr_info[1:]:
        pair = ' '.join(i.split()[2:4])
        info = i.split()[1]+'  '+' '.join(i.split()[4:6])
        opr[pair] = info
    count = 0
    for line in complete:
        pair = ' '.join(line.split()[2:4]) 
        if pair in opr:
            count += 1
            line = line.strip()+'  '+opr[pair]
        print(line.strip())
    print(count)
    
def add_column():
    with open('data/prokaryotic_gene_pairs/pairs_dataset2.txt') as file:
        data = file.readlines()
    # header = data[0].split()
    # new_data = [' '.join(header[:8])+' type '+' '.join(header[8:])]
    new_data = []
    for line in data:
        line = line.split()
        if len(line) == 8:
            line = ' '.join(line)+' Different_Operon'
        elif int(line[10]) - int(line[9]) == 1:
            line = ' '.join(line[:8])+' Adjacent '+' '.join(line[8:])
        elif int(line[10]) - int(line[9]) > 1:
            line = ' '.join(line[:8])+' Nonadjacent '+' '.join(line[8:])
        new_data.append(line)
    return new_data
        

def grep_operon(operon):
    dirs = glob.glob('./data/operon_datasets/*')
    dirs = [i for i in dirs if '_operon_' in i and 'rib' not in i]
    files = []
    for dirty in dirs:
        files += glob.glob(dirty+'/*')
    ops = []
    for filename in files:
        with open(filename) as file:
            for string in file:
                if operon in string:
                    ops.append(string+''+filename)
    count = 0
    for line in ops:
        if line.split()[0] == operon:
            fname = line.split()[-1]
            break
    new = [i for i in ops if i.split()[-1] == fname and i.split()[0] == operon]
    return(len(new))
    
def operon_counter():
    """Gets length of operons"""
    operons = []
    with open('data/prokaryotic_gene_pairs/pairs_w_dimers.txt') as file:
        data = file.readlines()
    for line in data:
        if len(line.split()) > 8:
            operons.append(line.split()[8])
    ops_dict = {op: grep_operon(op) for op in set(operons)}
    for line in data:
        if len(line.split()) > 8:
            line = line.strip()+' '+str(ops_dict[line.split()[8]])
        print(line.strip())
    return ops_dict


def old_data(filename):
    with open(filename) as file:
        data = file.readlines()
    op_dict = {}
    for line in data[1:]:
        line = line.split()
        op = line[1]
        p1 = line[4].strip('(,')
        p2 = line[5].strip(')')
        g1 = line[2]
        g2 = line[3]
        if op not in op_dict:
            op_dict[op] = {g1: p1, g2: p2}
        else:
            op_dict[op][g1] = p1
            op_dict[op][g2] = p2
    return op_dict
    
def new_data(filename, op_dict):
    with open(filename) as file:
        data = file.readlines()
    print('\t'.join(data[0].split()))
    for line in data[1:]:
        line = line.split()
        if len(line) > 8:
            g1 = line[2]
            g2 = line[3]
            op = line[8]
            new_pos = str((int(op_dict[op][g1]), int(op_dict[op][g2])))
            line = line[:9] + [new_pos] + [line[-1]]
        print('\t'.join(line))
            
od = old_data('data/int_pair2operon/combined_pairs_final.txt')
new_data('data/test.tmp', od)    
        
if __name__ =='__main__':
    od = old_data('data/int_pair2operon/combined_pairs_final.txt')
    new_data('data/test.tmp', od)    
    

    
