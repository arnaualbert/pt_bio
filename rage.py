import json
from   pprint   import pp
from   pathlib  import Path
from   Bio      import Entrez
from typing     import Union

import utils
from Bio.Entrez.Parser import DictionaryElement, ListElement
import utils

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Entrez.Parser import DictionaryElement, ListElement
from Bio.SeqIO.InsdcIO import GenBankIterator


Entrez.email = 'carde602@gmail.com'

CACHE_DIR = Path('/bio/practica_2/data')


# Main
# ---------------------------------------------------------------------
this_module: str = __name__
main_module: str = "__main__"

if this_module == main_module:

                #
    list_of_terms = ['Rabies lyssavirus isolate 18018LIB, complete genome','ebolavirus, complete genome','Marburg marburgvirus isolate MARV001, complete genome','Nipah virus, complete genome','Variola virus[ORGN]', 'TREPONEMA PALLIDUM TRIPLET ','Method of Immunization against the 4 serotypes of Dengue fever','Kyasanur forest disease virus isolate W6204 NS5 gene, partial cds']


    for i in list_of_terms:
        utils.request_search(db='nucleotide',term=i,retmax=1,xml_filename=f'/bio/practica_2/data/{i}.xml')
    
   # search_content: DictionaryElement = utils.read_xml('data/nucleotide-coronavirus-search.xml')
