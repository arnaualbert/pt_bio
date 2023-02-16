import json
from   pprint   import pp
from   pathlib  import Path
from   Bio      import Entrez
from typing     import Union
import re
import os
import glob
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

def count_files(path):
    file_list = glob.glob(path)
    num_of_files = len(file_list)
    return num_of_files
    


def get_xml_from_ncbi(terms):
    for i in terms:
        utils.request_search(db='nucleotide',term=i,retmax=1,xml_filename=f'/bio/practica_2/data/{i}.xml')


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

def rename_files_gb(list_of_gb):
        for file in list_of_gb:
            record_iter: GenBankIterator = SeqIO.parse(file, 'gb')
            for record in record_iter:
                if record.features:
                    for feature in record.features:
                        if feature.type == "source":
                            name_list = feature.qualifiers["organism"]
                            name = name_list[0]
                            reg_white_space = r'\s'
                            pat_white_space = re.compile(reg_white_space)
                            name = pat_white_space.sub('_', name)
                            reg_name = r'\/([^\/]+)$'
                            pat_name = re.compile(reg_name)
                            file_name = pat_name.search(file).group(1)
                            new_file_name = os.path.dirname(file) + '/' + name + os.path.splitext(file_name)[1]
                            os.rename(file, new_file_name)

# def rename_files(directory_path):
#     list_of_files = get_files_in_dir(directory_path)
#     list_of_parser_files = []
#     for file in list_of_files:
#         file_parsed:GenBankIterator = SeqIO.parse(file,'gb')
#         list_of_parser_files.append(file_parsed)
#     return list_of_parser_files

# Main
# ---------------------------------------------------------------------
this_module: str = __name__
main_module: str = "__main__"

if this_module == main_module:

    # List of term to search for the id
    list_of_terms = ['Rabies lyssavirus isolate 18018LIB, complete genome','ebolavirus, complete genome','Marburg marburgvirus isolate MARV001, complete genome','Nipah virus, complete genome','Variola virus[ORGN]', 'TREPONEMA PALLIDUM TRIPLET ','Method of Immunization against the 4 serotypes of Dengue fever','Kyasanur forest disease virus isolate W6204 NS5 gene, partial cds']
    
    data_to_count = '/bio/practica_2/data/*' 
    gb_files_to_count = 'practica_2/gb_files/*'
    data = '/bio/practica_2/data'
    gb_files = 'practica_2/gb_files'
    num_of_xml_files_before_search_in_ncbi = count_files(data_to_count)
    num_of_gb_files_before_search_in_ncbi = count_files(gb_files_to_count)
    # Get the xml from ncbi


    print(num_of_xml_files_before_search_in_ncbi)
    print(num_of_gb_files_before_search_in_ncbi)

    if num_of_xml_files_before_search_in_ncbi == 0:
        get_xml_from_ncbi(list_of_terms)
        print('excecuted: get_xml_from_ncbi')

    if num_of_xml_files_before_search_in_ncbi > 0:
        list_of_xml = get_files_in_dir(data)
        id_list = get_id_list(list_of_xml)
        print(id_list)
        get_gb_files(id_list)
        print('excecuted: get_gb_files')


    if num_of_gb_files_before_search_in_ncbi > 0:
        list_of_gb = get_files_in_dir(gb_files)
        rename_files_gb(list_of_gb)
        print('excecuted: rename_files_gb')
    

    # GET NAME OF THE ORGANISM
    # for rec in SeqIO.parse('/bio/practica_2/gb_files/EU293337.1.gb','gb'):
    #     if rec.features:
    #         for feature in rec.features:
    #             if feature.type == "source":
    #                 name_list = feature.qualifiers["organism"]
    #                 name = name_list[0]
    #                 print(name)
    #                 print(type(name))


    ################# GET CDS (MAKE A FUNCTION OG THIS) ###################################
    ##### https://stackoverflow.com/questions/23333123/extracting-cds-sequences-in-biopython
    # for rec in SeqIO.parse('/bio/practica_2/gb_files/EU293337.1.gb','gb'):
    #     if rec.features:
    #         for feature in rec.features:
    #             if feature.type == "CDS":
    #                 cds = feature.qualifiers["translation"]
    #                 print(cds)
    #                 print(type(cds))



    # features_list = []
    # for record in record_list:
    #     # print(record.features)
    #     # print()
    #     features_list.append(record.features)
    # print(features_list)

    # for features in features_list:
    #     if features.type == "CDS":
    #         print(features.location)
    #         print(features.qualifiers["protein_id"])
    #         print(features.location.extract(rec).seq)



    # parsed_files = rename_files(data)
    # record_list: list[SeqRecord] = list(parsed_files)
    # print(parsed_files)
    # for record in record_list:
    #     print(record)
    #     for rec in record:
    #         print(rec)