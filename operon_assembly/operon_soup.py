#!/usr/bin/env python3

"""Fairly certain this is redundant. save for legacy if you want."""

"""Program used to scrape DOOR2 database for operon information. 
14.09.2014"""

from bs4 import BeautifulSoup
import requests
import re

def open_page(url):
    req = requests.get(url)
    data = req.text
    return data
         
def parse_gen_info(data, ident=None):
    """Scrape the plasmid/chromosome operon file identifiers for any given organism strain."""
    soup = BeautifulSoup(data)
    re_obj = re.compile(r'([0-9]*) (NC_[0-9]*) (.*),')
    for item in soup.find_all('li')[8:]:
        info = re_obj.match(item.text)
        print("{:<8}".format(ident), end = '')
        print("{:<11}".format(info.group(1)), end = '')
        print("{:<16}".format(info.group(2)), end = '')
        print("{:<8}".format(info.group(3)), end = '\n')

def fetch_database_ids():
    """Loops through entirety of DOOR2s organisms strains in order to allow identification of operon filenames."""
    print("{:<8}{:<11}{:<16}{:<8}".format('GenId', 'Op_ID', 'NC_ID', 'Organism'))
    for x in range(1, 2080):
        try:
            url = 'http://csbl.bmb.uga.edu/DOOR/displayNC.php?id='+str(x)
            data = openPage(url)
            parseGenInfo(data, x)
        except:
            pass
        
if __name__ == "__main__":
    fetch_database_ids()