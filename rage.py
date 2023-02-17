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
from Bio.Seq        import Seq
from Bio.SeqFeature import SeqFeature
from Bio.Align      import PairwiseAligner, PairwiseAlignments, PairwiseAlignment, substitution_matrices
from Bio.Align.substitution_matrices import Array
import matplotlib.pyplot as plt

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

def get_alinaments_with_specific_virus(virus_to_compare,virus_list_to_compare):
    list_of_scores = []
    for virus in virus_list_to_compare:
        aligner = PairwiseAligner()
        score = aligner.score(virus_to_compare,virus)
        list_of_scores.append(score)
    return list_of_scores

def get_list_of_virus_to_compare(key_virus_good,complete_dict_of_virus_to_compare):
    print(len(list(complete_dict_of_virus_to_compare.values())))
    new_dict = complete_dict_of_virus_to_compare
    del new_dict[key_virus_good]
    list_virus_to_compare = new_dict
    print(len(list_virus_to_compare))
    return list_virus_to_compare

def get_alinaments_with_specific_virus_with_dict(virus_to_compare: str ,virus_dict_to_compare: dict):
    """
    Get the list of scores of the alignment of a virus to a virus
    input: virus_to_compare: str (virus to compare)
    output: list_of_scores: dict (list of scores)
    """
    list_of_socres = {}
    list_of_names_of_virus = list(virus_dict_to_compare.keys())
    list_of_cds_to_compare = list(virus_dict_to_compare.values())
    for index, name in enumerate(list_of_names_of_virus):
        aligner = PairwiseAligner()
        blosum62_matrix: Array = substitution_matrices.load('BLOSUM62')
        aligner.substitution_matrix = blosum62_matrix
        score = aligner.score(virus_to_compare,list_of_cds_to_compare[index])
        list_of_socres[name] = score
    return list_of_socres

def get_percentage(dict_of_scores):
    scores = list(dict_of_scores.values())
    max_value = max(scores)
    list_percentage = []
    for value in scores:
        if value == max_value:
            percentage = 100
            list_percentage.append(percentage)
        else:
            percentage = round((value/max_value)*100)
            list_percentage.append(percentage)
    return list_percentage

    

# Main
# ---------------------------------------------------------------------
this_module: str = __name__
main_module: str = "__main__"

if this_module == main_module:

    ### LIST OF TERMS ###

    list_of_terms = ['Homo sapiens insulin receptor (INSR) gene, complete cds','Macaca mulatta insulin receptor mRNA, partial cds, exons 9 - 12',
                     'Ovis aries partial mRNA for insulin receptor (ir gene)','Bos taurus partial mRNA for insulin receptor',
                     'Bubalus bubalis mRNA for insulin receptor, partial','Rattus norvegicus insulin receptor gene, partial exon 18, exon 19 and partial cds',
                     'Lepisma saccharina insulin receptor (InR) mRNA, partial cds','Bombyx mori insulin receptor (InR), mRNA',
                     'Oreochromis niloticus insulin receptor mRNA, partial cds']

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
        good_virus = 'Homo_sapiens'
        cds_virus = cds_dict['Homo_sapiens']
        socores = get_alinaments_with_specific_virus_with_dict(cds_virus,cds_dict)
        per = get_percentage(socores)
        results = pd.DataFrame()
        results['virus'] = list(cds_dict.keys())
        results['accession_number'] = accesion_numbers
        results['scores'] = list(socores.values()) 
        results['percentage'] = per
        results.to_csv('practica_2/results/infor.csv', index=False)
        plt.bar(results['virus'],results['percentage'])
        plt.savefig('practica_2/results/infor.png')
        print(results)
