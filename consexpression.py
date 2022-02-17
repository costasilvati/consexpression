#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Compreendendo o uso de argumentos - passo 2.
"""

import sys
import argparse


def main():
    parser = argparse.ArgumentParser(prog="Consexpression v1.1",
                                     description='Complete User Guide available: https://costasilvati.github.io/consexpression/',
                                     usage="")  # (1)
    parser.add_argument('--name', "-n",
                        required=True,
                        type=str,
                        help="Exepriment Name (text)")
    parser.add_argument('--samples', "-s",
                        required=False,
                        type=int,
                        help="Number of biological replics (integer value)")
    parser.add_argument('--groups', "-g",
                        required=True,
                        help="Number of groups to compare (integer value)")
    parser.add_argument('--groupname', "-gn",
                        required=True,
                        help="List with groups name (comma separetd list)")
    parser.add_argument('--reference', "-r",
                        required=True,
                        type=str,
                        help="Path to FASTA file, genome reference")
    parser.add_argument('--reads-path', "-rp",
                        required=True,
                        type=str,
                        help="Path to directory sequencing FASTQ file (reads)")
    parser.add_argument("--sub-dir", "-sd",
                        required=True,
                        type=list,
                        help="List for each group FASTQ directory")
    parser.add_argument('--threads', "-t",
                        type=int,
                        required=False,
                        default=2)
    parser.add_argument('--paired-end',"-pe",
                        type=str,
                        choices=['paired-end', 'single-end'],
                        required=True,
                        help="Single-end sequence use: single-end, Paired-end sequence, use: paired-end")
    parser.add_argument('--annotation', "-gtf",
                        type=str,
                        required=True,
                        help="Ptah to annotation GTF/GFF file")
    args = parser.parse_args()  # (3)

    print("Intervalo= {}".format(args.intervalo))  # (4)

    return 0


if __name__ == '__main__':
    sys.exit(main())