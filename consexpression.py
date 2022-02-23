import os.path
import sys
import argparse
import glob
from dao.experimentDao import ExperimentDao
from os.path import exists


def main():
    args = line_parameters()
    expDao = parameters_valid(args)
    if(expDao):
        exit(0)
    expDao._replic = get_reads_file(expDao._read_directory)
    # print (args)



def line_parameters():
    parser = argparse.ArgumentParser(prog="Consexpression v1.1\n",
        description='Complete User Guide available: https://costasilvati.github.io/consexpression/',
                                     usage="Developing...")  # (1)
    parser.add_argument('--name', "-n",required=True,type=str,
                        help="Exepriment Name (text)")
    parser.add_argument('--group_name', "-gn", nargs='+', type=str, required=True,
                        help="List with groups name, use a space separetd list.") #use to identify n of groups
    #parser.add_argument('--samples',"-s", nargs='+',type=int, required=True,
    #                help="Number of biological replics by group, use a integer value by each group(space separeted list).")
    parser.add_argument('--reference', "-r",type=str,
                        help="Path to FASTA file genome reference sequence.")
    parser.add_argument('--reads_path', "-rp",type=str,
                        help="Path sequencing FASTQ file (reads) to directory.")
    parser.add_argument("--sub_dir", "-sd", nargs='+',type=str,
                        help="List for each group FASTQ directory, need one directory by group.")
    parser.add_argument('--threads', "-t",type=int,required=False,default=2)
    parser.add_argument('--paired_end',"-pe", action='store_true',
                        help="Is a paired-end sequencing.")
    parser.add_argument('--annotation', "-gtf",type=str,required=True,
                        help="Ptah to annotation GTF/GFF file")
    parser.add_argument('--table_count', "-tb", type=str,
                        help="Path to file. Only you have mapped and count reads.")
    args = parser.parse_args()  # (3)
    print("\nIniciando Experimento {}".format(args.name))  # (4)
    print("\n")
    return args

def name_valid(name):
    """
    Verify the name: if is empty change to default name 'consexpression'
    :rtype: str
    :param name: name of experiment
    :return: str
    """
    if len(name) == 0:
        name = "consexpression"
        print("Experiment name is empty! The name was changed to consexpression")
    return name

def rep_valid(rep):
    """
    Verify if number of replicates technical and biological is valid (>= 1).
    int :param rep
    void :return
    """
    ok = False
    if rep >= 1:
        print("More then one!")
    else:
        print("1 replic or more (technique or biological)")
        print("number of replics in line 5 - 6")
        exit()

def parameters_valid(args):
    """
    Verify line parameters
    :param args: List of parameters
    :return: ExpreimentDao object
    """
    result_text = ""
    expDao = ExperimentDao()
    expDao._name = name_valid(args.name)

    if (args.table_count == None):
        print("Mapping, count and expression!")
        if (os.path.exists(args.reference) == False):
            result_text = result_text, "ERROR: \n - Reference file not found in path.", args.reference, "\n"
        else:
            expDao._reference = args.reference

        if (os.path.exists(args.annotation) == False):
            result_text = result_text, "ERROR: \n - Annotation file not found in path.", args.annotation, "\n"
        else:
            expDao._annotation_file = args.annotation
    else:
        print("Olnly expression...")

    expDao._group_name = args.group_name
    expDao._group_number = len(args.group_name)
    expDao._read_directory = args.reads_path
    expDao._replic = get_reads_file(expDao._read_directory)

    if(len(result_text) > 1 ):
        print(result_text)
        return None
    else:
        return expDao


def file_valid(self, path):
    """
    Verify if path is a file
    :param self:
    :param path: File path with extension
    :return: boolean
    """
    if os.path.isfile(path):
        print(" File: ", path, " found!")
        return True
    else:
        print(" ERROR: File: ", path, " NOT FOUND!")
        return False

def get_reads_file(listdir):
    """
    Get all fastq path of listdir
    :param listdir: List of directories by samples
    :return:
    """
    fastq_file = listdir
    length = len(listdir)
    i = 0
    # Iterating using while loop
    while i < length:
        path = listdir[i] + "*.fastq" #serach
        for file in glob.glob(path):
            fastq_file[listdir[i]].append(file)
            print(file)
        if len(fastq_file) == 0:
            print("ERROR: Not found files FASTQ in directory " + listdir[i])
    return fastq_file


if __name__ == '__main__':
    expDao = ExperimentDao()
    sys.exit(main())