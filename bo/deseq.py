# coding=utf-8

import rpy2.robjects as robjects
from bo.message import Message
from rpy2.rinterface import *
import warnings
warnings.filterwarnings("ignore", category=RRuntimeWarning)



class DESeq (object):
    """
    Run DESeq analysis
    """
    def __init__(self, count, group, repl, out):
        """
        Define the edgeR object
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
        self._logfc_column = 6
        self._pvalue_column = 7
        self._pvalue = 0.05
        self._logfc = 2

    def run_de(self, gene):
        de = 0
        try:
            lfc = float(gene[self._logfc_column])
            pv = float(gene[self._pvalue_column])
            if lfc >= self._logfc or lfc <= -self._logfc:
                if pv <= self._pvalue:
                    de = 1
        except ValueError:
            de = 0
        return de

    def run_deseq(self):
        """
        Execute default analysis with DESeq
        :return:
        """
        try:
            res = robjects.r('library("parallel")')
            res = robjects.r('library("stats4")')
            res = robjects.r('library("BiocGenerics")')
            res = robjects.r('library("Biobase")')
            res = robjects.r('library("locfit")')
            res = robjects.r('library(DESeq)')
            res = robjects.r('library("lattice")')
            ct = 'table <- read.csv("' + self._table_count + '",  row.names = 1, header = TRUE, stringsAsFactors=FALSE)'
            res = robjects.r(ct)
            res = robjects.r('m <- as.matrix(table)')
            grup = ""
            b_test = ""
            assert isinstance(self._replic, int)
            for ind in iter(self._groups_name):
                aux = "'" + ind + "', "
                b_test = aux + b_test
                grup = grup + aux * self._replic
            grup = grup[:(len(grup) - 2)]
            b_test = b_test[:len(b_test) - 2]
            res = robjects.r('condition = factor( c(' + grup + '))')
            res = robjects.r('cds <- newCountDataSet(m, condition)')
            res = robjects.r('cds <- estimateSizeFactors(cds)')
            command = ""
            if (self._replic == 1):
                command = 'cds <- estimateDispersions(cds, method="blind", fitType="local")' # fitType="local"
            else:
                command ='cds <- estimateDispersions(cds, fitType="local")' #fitType="local"

            res = robjects.r(command)
            cm = 'res <- nbinomTest(cds, ' + b_test + ')'
            res = robjects.r(cm)
            wr = 'write.table(res, file="' + self._output + '", sep = "\t", quote = FALSE)'
            res = robjects.r(wr)
        except RRuntimeError as rre:
            self._message.message_9("Error in DESeq execution: " + str(rre))
            #raise rre

        self._message.message_9("--- DESeq: is completed!")

# =============================== TESTES DA CLASSE ==================================
# inp = '/home/juliana/Documentos/Projeto_Juliana/Datasets/consexpression_replics_table_count.csv'
# gr = ["g1", "g2"]
# rp = 2
# out = '/home/juliana/Documentos/Projeto_Juliana/Datasets/consexpression_deseq.csv'
# t = DESeq(inp, gr, rp, out)
# t.run_deseq()