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
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Entrez.Parser import DictionaryElement, ListElement
from Bio.SeqIO.InsdcIO import GenBankIterator
from Bio import pairwise2


Entrez.email = 'carde602@gmail.com'

CACHE_DIR = Path('/bio/practica_2/data')
from pathlib import Path

def count_files(path):
    """
    Count the number of files in a directory
    Input: path: srt (path of the directory)
    Output: num_of_files: int (number of files in the directory)
    """
    file_list = glob.glob(path)
    num_of_files = len(file_list)
    return num_of_files
    
def get_xml_from_ncbi(terms):
    """
    Get the xml file from NCBI
    Input: terms: list (list of terms)
    Generate the xml files
    """
    for i in terms:
        utils.request_search(db='nucleotide',term=i,retmax=1,xml_filename=f'/bio/practica_2/data/{i}.xml')

def get_files_in_dir(directory_path):
    """
    Get the list of files in a directory
    Input: directory_path: str (path of the directory)
    Output: files: list (list of files in the directory)
    """
    directory_path = Path(directory_path)
    files = []
    for file in directory_path.glob('*'):
        if file.is_file():
            files.append(str(file))
    return files

def get_id_list(list_of_files):
    """
    Get the list of ids from a list of files
    Input: list_of_files: list (list of files)
    Output: id_list: list (list of ids)
    Read de xml and mekke it ain a dict where we extract the firs id of eachone
    """
    id_list = []
    for file in list_of_files:
        dict = utils.read_xml(file)
        id_for_list = dict['IdList'][0]
        id_list.append(id_for_list)
    return id_list

def get_gb_files(list_of_id):
    """
    Get the gb files from NCBI with a list of ids
    Input: list_of_id: list (list of ids)
    Output:  the genbank files
    """
    for id in list_of_id:
        utils.request_fetch('nucleotide',id,'gb',f'/bio/practica_2/gb_files/{id}.gb')

def rename_files_gb(list_of_gb):
        """
        Rename the genbank files with regex 
        Input: list_of_gb: list (list of genbank files)
        Rename the genbank files
        """
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
                            if file != new_file_name:
                                os.rename(file, new_file_name)
                            else:
                                print(f'The file {new_file_name} have the good name')

def get_dict_of_cds(gb_files):
        """
        Get the dict of cds from a list of genbank files
        Input: gb_files: list (list of genbank files)
        Output: cds_dict: dict (dict of cds)
        """
        list_of_gb = get_files_in_dir(gb_files)
        print(list_of_gb)
        cds_list = []
        cds_dict = {}
        for file in list_of_gb:
            record_iter: GenBankIterator = SeqIO.parse(file, 'gb')
            for record in record_iter:
                if record.features:
                    list_of_cds_seq = []
                    for feature in record.features:
                        if feature.type == "source":
                            name_list = feature.qualifiers["organism"]
                            name = name_list[0]
                            reg_white_space = r'\s'
                            pat_white_space = re.compile(reg_white_space)
                            name = pat_white_space.sub('_', name)
                        if feature.type == "CDS":
                            cds = feature.qualifiers["translation"]
                            list_of_cds_seq.append(cds)
                    print(f'{name} have: {len(list_of_cds_seq)} CDS')
                    cds_list.append(list_of_cds_seq[0])
                    list_elem_one = list_of_cds_seq[0]
                    cds = list_elem_one[0]
                    cds_dict[name] = cds
        return cds_dict

def get_accession_number(gb_files):
    """
    get the accession number from a list of genbank files
    input: gb_files: list (list of genbank files)
    output: accession_numbers_list: list (list of accession numbers)
    """
    list_of_gb = get_files_in_dir(gb_files)
    accesion_numbers_list = []
    for file in list_of_gb:
        record_iter: GenBankIterator = SeqIO.parse(file, 'gb')
        for record in record_iter:
            if record.name:
                accesion_numbers_list.append(record.name)
    return accesion_numbers_list

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

    ### LIST OF TERMS ###

    list_of_terms = ['Rabies lyssavirus isolate 18018LIB, complete genome','ebolavirus, complete genome','Marburg marburgvirus isolate MARV001, complete genome','Nipah virus, complete genome','Variola virus[ORGN]', 'TREPONEMA PALLIDUM TRIPLET ','Kyasanur forest disease virus isolate W6204 NS5 gene, partial cds','Hepatitis B virus isolate I172, complete genome','Human rhinovirus B strain KR2629 polyprotein gene, partial cds','Adenovirus type 2, complete genome']

    ### PATH OF THE FOLDER THAT CONTAINS THE XML FILES AND THE GB FILES ###

    data_to_count = '/bio/practica_2/data/*' 
    gb_files_to_count = 'practica_2/gb_files/*'
    data = '/bio/practica_2/data'
    gb_files = 'practica_2/gb_files'

    ### COUNT THE NUMBER OF XML FILES AND GB FILES ###

    num_of_xml_files_before_search_in_ncbi = count_files(data_to_count)
    num_of_gb_files_before_search_in_ncbi = count_files(gb_files_to_count)

    ### INFORMS THE USER THE NNUMBER OF XML FILES AND GB FILES ###

    print(num_of_xml_files_before_search_in_ncbi)
    print(num_of_gb_files_before_search_in_ncbi)

    ### IF THE USER DON'T HAVE THE XML FILES THAT HE WANTS THE FUNCTION IS EXECUTED ###

    if num_of_xml_files_before_search_in_ncbi == 0:
        get_xml_from_ncbi(list_of_terms)
        print('excecuted: get_xml_from_ncbi')

    ### IF THE USER DON'T HAVE THE GB FILES THAT HE WANTS THE THE FUNCTIONS ARE EXECUTED ###

    if num_of_gb_files_before_search_in_ncbi == 0:
        list_of_xml = get_files_in_dir(data)
        id_list = get_id_list(list_of_xml)
        print(id_list)
        get_gb_files(id_list)
        print('excecuted: get_gb_files')
        list_of_gb = get_files_in_dir(gb_files)
        rename_files_gb(list_of_gb)
        print('excecuted: rename_files_gb')

    ### THE USER GETS A PANDAS DATAFRAME WITH THE INFORMATION OF THE GENBANK FILES ###

    if num_of_gb_files_before_search_in_ncbi != 0:
        cds_dict = get_dict_of_cds(gb_files)
        accesion_numbers = get_accession_number(gb_files)
        results = pd.DataFrame()
        results['virus'] = list(cds_dict.keys())
        results['accession_number'] = accesion_numbers 
        print(results)



    # print(cds_dict.keys)
    # for key in cds_dict.keys():
    #     print(key)
    # list_of_gb = get_files_in_dir(gb_files)
    # for file in list_of_gb:
    #     record_iter: GenBankIterator = SeqIO.parse(file, 'gb')
    #     for record in record_iter:
    #         if record.name:
    #             print(record.name)






#############################old get dic of cds#############################
        # print(cds_dict)
        # list_of_gb = get_files_in_dir(gb_files)
        # cds_list = []
        # cds_dict = {}
        # for file in list_of_gb:
        #     record_iter: GenBankIterator = SeqIO.parse(file, 'gb')
        #     for record in record_iter:
        #         if record.features:
        #             list_of_cds_seq = []
        #             for feature in record.features:
        #                 if feature.type == "source":
        #                     name_list = feature.qualifiers["organism"]
        #                     name = name_list[0]
        #                     reg_white_space = r'\s'
        #                     pat_white_space = re.compile(reg_white_space)
        #                     name = pat_white_space.sub('_', name)
        #                 if feature.type == "CDS":
        #                         # print(feature.location)
        #                         # print(feature.qualifiers["protein_id"])
        #                         # print(feature.location.extract(record).seq)
        #                         # list_of_cds_seq = []
        #                     cds = feature.qualifiers["translation"]
        #                     list_of_cds_seq.append(cds)
        #                     # cds_list.append(cds)
        #             print(f'{name} have: {len(list_of_cds_seq)} CDS')
        #             cds_list.append(list_of_cds_seq[0])
        #             list_elem_one = list_of_cds_seq[0]
        #             cds = list_elem_one[0]
        #             cds_dict[name] = cds
        #             return cds_dict
        # print(cds_list)
        # print(len(cds_list))
        # print(len(cds_dict))






    # if num_of_gb_files_before_search_in_ncbi > 0:
        # list_of_gb = get_files_in_dir(gb_files)
        # rename_files_gb(list_of_gb)
        # print('excecuted: rename_files_gb')

    # list_of_gb = get_files_in_dir(gb_files)
    # for file in list_of_gb:
    #     cds_list = []
    #     record_iter: GenBankIterator = SeqIO.parse(file, 'gb')
    #     for record in record_iter:
    #         if record.features:
    #             for feature in record.features:
    #                 if feature.type == "CDS":
    #                     cds = feature.qualifiers["translation"][0]
    #                     cds_list.append(cds)
    # print(cds_list)

    

    



 
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