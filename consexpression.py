import os.path
import sys
import argparse
import glob
from dao.experimentDao import ExperimentDao
from os.path import exists


def main():
    args = line_parameters()
    exp_dao = parameters_valid(args)
    if(exp_dao):
        exit(0)
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
    parser.add_argument('--reads_path', "-rp",nargs='+', type=str,
                        help="Path sequencing FASTQ file (reads) to directory (one for each treatment)")
    # parser.add_argument("--sub_dir", "-sd", nargs='+',type=str,
    #                     help="List for each group FASTQ directory, need one directory by group.")
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
        print("More then one replic!")
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
    exp_dao_tmp = ExperimentDao()
    exp_dao_tmp._name = name_valid(args.name)
    exp_dao_tmp._group_name = args.group_name
    try:
        if (args.table_count == None):
            print("Mapping, count and expression!")
            if (os.path.exists(args.reference) == False):
                result_text = result_text, "ERROR: \n - Reference file not found in path.", args.reference, "\n"
            else:
                exp_dao_tmp._reference = args.reference

            if (os.path.exists(args.annotation) == False):
                result_text = result_text, "ERROR: \n - Annotation file not found in path.", args.annotation, "\n"
            else:
                exp_dao_tmp._annotation_file = args.annotation

            if (os.path.exists(args.reads_path)):
                exp_dao_tmp._read_directory = args.reads_path
                exp_dao_tmp._replic = get_reads_file(expDao._read_directory)

                if(exp_dao_tmp._replic == None):
                    result_text = result_text, "FASTQ reads not found in ", exp_dao_tmp._read_directory, "\n"
                else:
                    relation_file_group(exp_dao_tmp._replic, exp_dao_tmp._group_name)
            else:
                result_text = result_text, "Path ",args.reads_path," to reads directory NOT FOUND. \n"
        else:
            print("Olnly expression...")
    except TypeError:
        print("ERROR: Some parameter have a wrong type!")
    except:
        print("Unspected exception!")


    if(len(result_text) > 2 ):
        print(result_text)
        return None
    else:
        return exp_dao_tmp


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

def get_reads_file(path_dir):
    """
    Get all fastq path of listdir
    :param listdir: List of directories by samples
    :return:
    """
    fastq_files = None
    path = path_dir + "*.fastq" #serach
    for file in glob.glob(path):
        fastq_files.append(file)
        print("File: ", file , "found!")
    if (fastq_files != None):
        pass
    else:
        print("ERROR: Not found files FASTQ in directory " + path_dir)

    return fastq_files


def relation_file_group (fastq_files, group_names):
    print("--- Organize FASTQ files by Treatment ----")
    print("Group numbers:")
    print(group_names)
# STOP - test fastq found


if __name__ == '__main__':
    expDao = ExperimentDao()
    sys.exit(main())