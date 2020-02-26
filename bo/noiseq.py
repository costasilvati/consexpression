# coding=utf-8

import rpy2.robjects as robjects
from rpy2.rinterface_lib.embedded import RRuntimeError

from bo.message import Message
from rpy2.rinterface import *
import warnings
warnings.filterwarnings("ignore", category=RRuntimeWarning)


class Noiseq(object):

    def __init__(self, count, group, repl, out):
        """
        Define the NOISeq object
        :param count:
        :param group:
        :param repl:
        :param out:
        """
        self._table_count = count
        self._groups_name = group
        self._replic = repl
        self._output = out
        self._message = Message()
        self._likelihood_column = len(group) + 3
        self._likelihood = 0.95

    def run_de(self, gene):
        de = 0
        try:
            like = float(gene[self._likelihood_column])
            if like >= self._likelihood:
                de = 1
        except ValueError:
            de = 0
        return de


    def run_noiseq(self):
        """
        Execute default analysis with NOISeq
        :return:
        """
        try:
            res = robjects.r('library("parallel")')
            res = robjects.r('library("splines")')
            res = robjects.r('library("Matrix")')
            res = robjects.r('library("BiocGenerics")')
            res = robjects.r('library("Biobase")')
            res = robjects.r('library("NOISeq")')
            ct = 'table <- read.csv("' + self._table_count + '",  row.names = 1, header = TRUE, stringsAsFactors=FALSE)'
            res = robjects.r(ct)
            res = robjects.r('table <- as.matrix(table)')
            ts = ""
            run = ""
            tsrun = ""
            count_run = 1
            assert isinstance(self._replic, int)
            for ind in iter(self._groups_name):
                aux = "'" + ind + "', "
                ts = ts + aux * self._replic
                while (count_run <= self._replic):
                    tsrun = tsrun + "'" + ind + str(count_run) + "', "
                    run = run + "'" + "R" + str(count_run) + "', "
                    count_run += 1
                count_run = 1
            ts = ts[:(len(ts) - 2)]
            tsrun = tsrun[:(len(tsrun) - 2)]
            run = run[:(len(run) - 2)]
            res = robjects.r('myfactors = data.frame(Tissue=c('+ ts +'), TissueRun=c(' + tsrun + '), Run=c(' + run + '))')
            res = robjects.r('mydata <- readData(data = table, factors = myfactors)')
            res = robjects.r('mynoiseq = noiseq(mydata, k = 0.5, factor = "Tissue", lc = 1, replicates = "technical")')
            res = robjects.r('results <- head(mynoiseq@results)')
            res = robjects.r('write.csv(results, file="' + self._output + '", sep = "\t", quote = FALSE)')
            self._message.message_9("--- NOISeq: is completed!")
        except RRuntimeError as rre:
            self._message.message_9("Error in NOISeq execution: " + str(rre))
            raise rre
#========================= TESTE da CLASSE==============
# inp = 'UHR_vs_Brain_gencode_TopHat_NOISeq.csv'
# inp = 'consexpression_NOISeq.csv'
# grp = "g1", "g2"
# rep = 1
# out = 'consexpression_NOISeq_out.csv'
# b = Noiseq(inp, grp, rep, out)
# read_bay = open(inp, 'r')
# c_b = 0
# for line in iter(read_bay):
#     #print('--' + line)
#     if c_b > 0:
#        gene = line.split(",")
#        print(gene[0])
#        v = b.run_de(gene)
#        print('--> '+ str(v))
#     c_b += 1