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
from pathlib import Path

def get_files_in_dir(directory_path):
    directory_path = Path(directory_path)
    files = []
    for file in directory_path.glob('*'):
        if file.is_file():
            files.append(str(file))
    return files

def get_id_list(list_of_files):
    id_list = []
    for file in list_of_files:
        dict = utils.read_xml(file)
        id_for_list = dict['IdList'][0]
        id_list.append(id_for_list)
    return id_list

def get_gb_files(list_of_id):
    for id in list_of_id:
        utils.request_fetch('nucleotide',id,'gb',f'/bio/practica_2/gb_files/{id}.gb')

def rename_files(directory_path):
    list_of_files = get_files_in_dir(directory_path)
    list_of_parser_files = []
    for file in list_of_files:
        file_parsed:GenBankIterator = SeqIO.parse(file,'gb')
        list_of_parser_files.append(file_parsed)
    return list_of_parser_files

# Main
# ---------------------------------------------------------------------
this_module: str = __name__
main_module: str = "__main__"

if this_module == main_module:

    # List of term to search for the id
    ###############     list_of_terms = ['Rabies lyssavirus isolate 18018LIB, complete genome','ebolavirus, complete genome','Marburg marburgvirus isolate MARV001, complete genome','Nipah virus, complete genome','Variola virus[ORGN]', 'TREPONEMA PALLIDUM TRIPLET ','Method of Immunization against the 4 serotypes of Dengue fever','Kyasanur forest disease virus isolate W6204 NS5 gene, partial cds']
    
    # Folder to save the id
    data = '/bio/practica_2/data' 

    # List of files 
    ############## list_of_files = get_files_in_dir(data)

    # List of id to search for in the ncbi
    ###############        id_list = get_id_list(list_of_files)

    # Get the gb files from the list of id
    ############    get_gb_files(id_list)


    record_iter: GenBankIterator = SeqIO.parse('/bio/practica_2/gb_files/EU293337.1.gb', 'gb')
    record_list: list[SeqRecord] = list(record_iter)

    for record in record_list:
        print(record)
        print()

    # parsed_files = rename_files(data)
    # record_list: list[SeqRecord] = list(parsed_files)
    # print(parsed_files)
    # for record in record_list:
    #     print(record)
    #     for rec in record:
    #         print(rec)