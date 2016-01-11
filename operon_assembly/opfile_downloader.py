#!/usr/bin/env python3

# just print organisms and their occurences. 15.09.2014.

import re
import collections
from bs4 import BeautifulSoup
import requests
import urllib.request


class GetStrains(object):
    """Gets the all the unique organisms (structures) and finds strain info.

    Arguments: A filename in the *_pairs.txt format.
    """

    def __init__(self, filename):
        with open(filename) as file:
            self.data = file.readlines()

    def get_strucs(self):
        strucs = []
        for line in self.data:
            strucs.append(line.split()[0].split('_')[0])
            strucs.append(line.split()[1].split('_')[0])
        return set(strucs)

    def parse_strains(self):
        """Returns a list of structures from the most common organisms."""
        line_pattern = re.compile(r'([a-z0-9]*)\s([a-z_]*)')
        org_list = []
        org_dict = {}
        structure_list = []
        for line in self.data:
            hit = line_pattern.match(line)
            structure = hit.group(1)
            organism = hit.group(2)
            org_list.append(organism)
            if structure not in structure_list:
                structure_list.append(structure)
            if organism not in org_dict:
                org_dict[organism] = [structure]
            else:
                org_dict[organism].append(structure)
        return structure_list

    def query_pdb(self,structure_list, printed=False):
        """Query PDB with list of structures to acquire specific strain info.
        Shoddy way of doing this. No need to scrape, just use API!
        """
        url = 'http://www.rcsb.org/pdb/rest/describeMol?structureId='
        taxonomy_dict = {}
        for item in structure_list:
            pdb = str(item)
            request = requests.get(url+pdb)
            data = request.text
            soup = BeautifulSoup(data)
            # print(item+' done')
            # Parse the specific strain names and id for each structure
            for info in soup.find_all('taxonomy'):
                tax_info = re.search(r'id="(\d*)" name="(.*)">', str(info))
                tax_id = tax_info.group(1)
                name = tax_info.group(2)
                if name not in taxonomy_dict:
                    taxonomy_dict[name] = tax_id
        if printed == True:
            for item in taxonomy_dict.items():
                print(item[0], end='    ')
                print(item[1])
        return taxonomy_dict

    def id_grabber(self, taxonomy_dict):
        """Returns a set of ids corresponding to operon files in DOOR2."""
        with open('data/int_pair2operon/door_gen_ids.txt') as door_gen_file:
            door_data = door_gen_file.readlines()
        door_pattern = re.compile(r'(\d*)\s*(\d*)\s*([A-Z_0-9]*)\s*(.*)')
        opr_id_list = []
        for line in door_data:
            strain_hit = door_pattern.match(line).group(4)
            opfile_id = door_pattern.match(line).group(2)
            for name in taxonomy_dict:
                try:
                    # Matches strains from tax_inf_file with those from
                    # door_gen_ids
                    succ_match = re.search(name, strain_hit)
                    if succ_match != None:
                        opr_id_list.append((opfile_id, strain_hit,
                                            taxonomy_dict[name]))
                except:
                    continue
        return set(opr_id_list)


def file_downloader(opr_id_set, dl_location):
    """Does the grunt work of mass downloading operon files from list of ids.
    """
    base_url = 'http://csbl.bmb.uga.edu/DOOR/downloadNCoperon.php?NC_id='
    for item in opr_id_set:
        url = base_url+item
        full_file_loc = dl_location+item+'.opr'
        try:
            urllib.request.urlretrieve(url, full_file_loc)
        except:
            print('could not download: '+url)
            continue
    print('done')

if __name__ == "__main__":
    pass