#!/usr/bin/env python

# HTSeq_QA.py
#
# (c) Simon Anders, European Molecular Biology Laboratory, 2010
# released under GNU General Public License

import sys
import os.path
import argparse
from itertools import islice
import numpy as np
import HTSeq

try:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import Normalize
except ImportError:
    sys.stderr.write("htseq-qa needs 'matplotlib >= 1.5'")
    raise


def get_read_length(readfile, isAlnmntFile):
    readlen = 0
    if isAlnmntFile:
        reads = (a.read for a in readfile)
    else:
        reads = readfile
    for r in islice(reads, 10000):
        if len(r) > readlen:
            readlen = len(r)

    return readlen


def compute_quality(
        readfilename,
        file_type,
        nosplit,
        readlen,
        max_qual,
        gamma,
        primary_only=False,
        max_records=-1,
        ):

    if file_type in ("sam", "bam"):
        readfile = HTSeq.BAM_Reader(readfilename)
        isAlnmntFile = True
    elif file_type == "solexa-export":
        readfile = HTSeq.SolexaExportReader(readfilename)
        isAlnmntFile = True
    elif file_type == "fastq":
        readfile = HTSeq.FastqReader(readfilename)
        isAlnmntFile = False
    elif file_type == "solexa-fastq":
        readfile = HTSeq.FastqReader(readfilename, "solexa")
        isAlnmntFile = False
    else:
        raise ValueError('File format not recognized: {:}'.format(file_type))

    twoColumns = isAlnmntFile and (not nosplit)

    if readlen is None:
        readlen = get_read_length(readfile, isAlnmntFile)

    # Initialize count arrays
    base_arr_U = np.zeros((readlen, 5), np.int64)
    qual_arr_U = np.zeros((readlen, max_qual+1), np.int64)
    if twoColumns:
        base_arr_A = np.zeros((readlen, 5), np.int64)
        qual_arr_A = np.zeros((readlen, max_qual+1), np.int64)

    # Main counting loop
    i = 0
    try:
        for a in readfile:
            if isAlnmntFile:
                r = a.read
            else:
                r = a

            # Exclude non-primary alignments if requested
            if isAlnmntFile and primary_only:
                if a.aligned and a.not_primary_alignment:
                    continue

            if twoColumns and isAlnmntFile and a.aligned:
                r.add_bases_to_count_array(base_arr_A)
                r.add_qual_to_count_array(qual_arr_A)
            else:
                r.add_bases_to_count_array(base_arr_U)
                r.add_qual_to_count_array(qual_arr_U)

            i += 1

            if i == max_records:
                break

            if (i % 200000) == 0:
                if (not isAlnmntFile) or primary_only:
                    print(i, "reads processed")
                else:
                    print(i, "alignments processed")

    except:
        sys.stderr.write("Error occured in: %s\n" %
                         readfile.get_line_number_string())
        raise

    if (not isAlnmntFile) or primary_only:
        print(i, "reads processed")
    else:
        print(i, "alignments processed")

    # Normalize result
    def norm_by_pos(arr):
        arr = np.array(arr, np.float64)
        arr_n = (arr.T / arr.sum(1)).T
        arr_n[arr == 0] = 0
        return arr_n

    def norm_by_start(arr):
        arr = np.array(arr, np.float64)
        arr_n = (arr.T / arr.sum(1)[0]).T
        arr_n[arr == 0] = 0
        return arr_n

    result = {
        'isAlnmntFile': isAlnmntFile,
        'readlen': readlen,
        'twoColumns': twoColumns,
        'base_arr_U_n': norm_by_pos(base_arr_U),
        'qual_arr_U_n': norm_by_start(qual_arr_U),
        'nreads_U': base_arr_U[0, :].sum(),
        }

    if twoColumns:
        result['base_arr_A_n'] = norm_by_pos(base_arr_A)
        result['qual_arr_A_n'] = norm_by_start(qual_arr_A)
        result['nreads_A'] = base_arr_A[0, :].sum()

    return result


def plot(
        result,
        readfilename,
        outfile,
        max_qual,
        gamma,
        primary_only=False,
        ):

    def plot_bases(arr, ax):
        xg = np.arange(readlen)
        ax.plot(xg, arr[:, 0], marker='.', color='red')
        ax.plot(xg, arr[:, 1], marker='.', color='darkgreen')
        ax.plot(xg, arr[:, 2], marker='.', color='lightgreen')
        ax.plot(xg, arr[:, 3], marker='.', color='orange')
        ax.plot(xg, arr[:, 4], marker='.', color='grey')
        ax.set_xlim(0, readlen-1)
        ax.set_ylim(0, 1)
        ax.text(readlen*.70, .9, "A", color="red")
        ax.text(readlen*.75, .9, "C", color="darkgreen")
        ax.text(readlen*.80, .9, "G", color="lightgreen")
        ax.text(readlen*.85, .9, "T", color="orange")
        ax.text(readlen*.90, .9, "N", color="grey")

    if outfile is None:
        outfilename = os.path.basename(readfilename) + ".pdf"
    else:
        outfilename = outfile

    isAlnmntFile = result['isAlnmntFile']
    readlen = result['readlen']
    twoColumns = result['twoColumns']

    base_arr_U_n = result['base_arr_U_n']
    qual_arr_U_n = result['qual_arr_U_n']
    nreads_U = result['nreads_U']

    if twoColumns:
        base_arr_A_n = result['base_arr_A_n']
        qual_arr_A_n = result['qual_arr_A_n']
        nreads_A = result['nreads_A']

    cur_backend = matplotlib.get_backend()

    try:
        matplotlib.use('PDF')

        fig = plt.figure()
        fig.subplots_adjust(top=.85)
        fig.suptitle(os.path.basename(readfilename), fontweight='bold')

        if twoColumns:

            ax = fig.add_subplot(221)
            plot_bases(base_arr_U_n, ax)
            ax.set_ylabel("proportion of base")
            ax.set_title(
                    "non-aligned reads\n{:.0%} ({:.4f} million)".format(
                    1.0 * nreads_U / (nreads_U+nreads_A),
                    1.0 * nreads_U / 1e6,
                    ))

            ax2 = fig.add_subplot(222)
            plot_bases(base_arr_A_n, ax2)
            ax2.set_title(
                    "{:}\n{:.0%} ({:.4f} million)".format(
                        'aligned reads' if primary_only else 'alignments',
                        1.0 * nreads_A / (nreads_U+nreads_A),
                        1.0 * nreads_A / 1e6,
                    ))

            ax3 = fig.add_subplot(223)
            ax3.pcolor(
                    qual_arr_U_n.T ** gamma,
                    cmap=plt.cm.Greens,
                    norm=Normalize(0, 1))
            ax3.set_xlim(0, readlen-1)
            ax3.set_ylim(0, max_qual+1)
            ax3.set_xlabel("position in read")
            ax3.set_ylabel("base-call quality score")

            ax4 = fig.add_subplot(224)
            ax4.pcolor(
                    qual_arr_A_n.T ** gamma,
                    cmap=plt.cm.Greens,
                    norm=Normalize(0, 1))
            ax4.set_xlim(0, readlen-1)
            ax4.set_ylim(0, max_qual+1)
            ax4.set_xlabel("position in read")

        else:

            ax = fig.add_subplot(211)
            plot_bases(base_arr_U_n, ax)
            ax.set_ylabel("proportion of base")
            ax.set_title("{:.3f} million {:}".format(
                1.0 * nreads_U / 1e6,
                'reads' if (not isAlnmntFile) or primary_only else 'alignments',
                ))

            ax2 = fig.add_subplot(212)
            ax2.pcolor(
                    qual_arr_U_n.T ** gamma,
                    cmap=plt.cm.Greens,
                    norm=Normalize(0, 1))
            ax2.set_xlim(0, readlen-1)
            ax2.set_ylim(0, max_qual+1)
            ax2.set_xlabel("position in read")
            ax2.set_ylabel("base-call quality score")

        fig.savefig(outfilename)

    finally:
        matplotlib.use(cur_backend)


def main():

    # **** Parse command line ****
    pa = argparse.ArgumentParser(
        description=
        "This script take a file with high-throughput sequencing reads " +
        "(supported formats: SAM, Solexa _export.txt, FASTQ, Solexa " +
        "_sequence.txt) and performs a simply quality assessment by " +
        "producing plots showing the distribution of called bases and " +
        "base-call quality scores by position within the reads. The " +
        "plots are output as a PDF file.",
        )
    pa.add_argument(
        'readfilename',
        help='The file to count reads in (SAM/BAM or Fastq)',
        )
    pa.add_argument(
        "-t", "--type", type=str, dest="type",
        choices=("sam", "bam", "solexa-export", "fastq", "solexa-fastq"),
        default="sam", help="type of read_file (one of: sam [default], bam, " +
        "solexa-export, fastq, solexa-fastq)")
    pa.add_argument(
        "-o", "--outfile", type=str, dest="outfile",
        help="output filename (default is <read_file>.pdf)")
    pa.add_argument(
        "-r", "--readlength", type=int, dest="readlen",
        help="the maximum read length (when not specified, the script guesses from the file")
    pa.add_argument(
        "-g", "--gamma", type=float, dest="gamma",
        default=0.3,
        help="the gamma factor for the contrast adjustment of the quality score plot")
    pa.add_argument(
        "-n", "--nosplit", action="store_true", dest="nosplit",
        help="do not split reads in unaligned and aligned ones")
    pa.add_argument(
        "-m", "--maxqual", type=int, dest="maxqual", default=41,
        help="the maximum quality score that appears in the data (default: 41)")
    pa.add_argument(
        '--primary-only', action='store_true',
        help="For SAM/BAM input files, ignore alignments that are not primary. " +
        "This only affects 'multimapper' reads that align to several regions " +
        "in the genome. By choosing this option, each read will only count as " +
        "one; without this option, each of its alignments counts as one."
    )
    pa.add_argument(
        '--max-records', type=int, default=-1, dest='max_records',
        help="Limit the analysis to the first N reads/alignments."
        )

    args = pa.parse_args()

    result = compute_quality(
        args.readfilename,
        args.type,
        args.nosplit,
        args.readlen,
        args.maxqual,
        args.gamma,
        args.primary_only,
        args.max_records,
        )

    plot(
        result,
        args.readfilename,
        args.outfile,
        args.maxqual,
        args.gamma,
        args.primary_only,
        )


if __name__ == "__main__":
    main()
