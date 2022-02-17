import sys
import argparse
from dao.experimentDao import ExperimentDao


def main():
    args = line_parameters()

    expDao = ExperimentDao()
    expDao._name = name_valid(args.name)
    # --- Validar num. replicatas - STOP
    expDao._replic = args.samples

    if (args.table_count == None):
        print ("Mapping, count and expression!")
    else:
        print("Olnly expression...")
    # print (args)



def line_parameters():
    parser = argparse.ArgumentParser(prog="Consexpression v1.1\n",
        description='Complete User Guide available: https://costasilvati.github.io/consexpression/',
                                     usage="Developing...")  # (1)
    parser.add_argument('--name', "-n",required=True,type=str,
                        help="Exepriment Name (text)")
    parser.add_argument('--groupname', "-gn", nargs='+', type=str, required=True,
                        help="List with groups name, use a space separetd list.") #use to identify n of groups
    parser.add_argument('--samples',"-s", nargs='+',type=int, required=True,
                    help="Number of biological replics by group, use a integer value by each group(space separeted list).")
    parser.add_argument('--reference', "-r",type=str,
                        help="Path to FASTA file genome reference sequence.")
    parser.add_argument('--reads-path', "-rp",type=str,
                        help="Path sequencing FASTQ file (reads) to directory.")
    parser.add_argument("--sub-dir", "-sd", nargs='+',type=str,
                        help="List for each group FASTQ directory, need one directory by group.")
    parser.add_argument('--threads', "-t",type=int,required=False,default=2)
    parser.add_argument('--paired-end',"-pe", action='store_true',
                        help="Is a paired-end sequencing.")
    parser.add_argument('--annotation', "-gtf",type=str,required=True,
                        help="Ptah to annotation GTF/GFF file")
    parser.add_argument('--table-count', "-tb", type=str,
                        help="Path to file. Only you have mapped and count reads.")
    args = parser.parse_args()  # (3)
    print("\nIniciando Experimento {}".format(args.name))  # (4)
    print("\n")
    return args

def name_valid(name):
    """
    Verify the name: if is empty change to default name
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


if __name__ == '__main__':
    expDao = ExperimentDao()
    sys.exit(main())