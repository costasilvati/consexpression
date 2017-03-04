# coding=utf-8

from rpy2.robjects import r
import rpy2.robjects as robjects
from bo.message import Message
from rpy2.rinterface import *
import warnings
warnings.filterwarnings("ignore", category=RRuntimeWarning)

class EdgeR(object):
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
        self._column_result = [3,4]
        self._min_result = []
        self._message = Message()
        self._logfc_colum = 1
        self._pvalue_colum = 3
        self._pvalue = 0.05
        self._logfc = 2


    def run_de(self, gene):
        de = 0
        lfc = float(gene[self._logfc_colum])
        pv = float(gene[self._pvalue_colum])
        if lfc >= self._logfc or lfc <= -self._logfc:
            if pv >= self._pvalue:
                de = 1
        return de

    def run_edger(self):
        """
        Execute default analysis with edegeR
        :return:
        """
        try:
            finish_message = ""
            res = robjects.r('library("limma")')
            res = robjects.r('library("edgeR")')
            ct = 'table <- read.csv("' \
                 + self._table_count + '",  row.names = 1, header = TRUE, stringsAsFactors=FALSE, sep = "' + "," + '")'
            res = robjects.r(ct)
            res = robjects.r('m <- as.matrix(table)')
            grup = ""
            assert isinstance(self._replic, int)
            for ind in iter(self._groups_name):
                aux = "'" + ind + "', "
                grup = grup + aux * self._replic
            grup = grup[:(len(grup) - 2)]
            grup = 'group <- c(' + grup + ')'
            res = robjects.r(grup)
            res = robjects.r('y.dge <- DGEList(counts = m, group = group)')
            if (self._replic < 1):
                self._message.message_4(" Replicates not found by edgeR. EdgeR should be executed manual form.")
            elif (self._replic == 1):
                # edgeR manual based solution for without replicates
                res = robjects.r('bcv <- 0.2')
                res = robjects.r('y.et <- exactTest(y.dge, dispersion = bcv^2)')
                res = robjects.r('y.tp <- topTags(y.et, n = 100000)')
                res = robjects.r('y.pvalues <- y.et$table$PValue')
                wr = 'write.table(y.tp$table, "' + self._output + '", sep = "\t", quote = FALSE)'
                res = robjects.r(wr)
                finish_message = "--- edgeR without replicates is completed!"
            else:
                r('y.dge <- calcNormFactors(y.dge)')
                r('y.dge <- estimateDisp(y.dge)')
                r('y.dge <- estimateCommonDisp(y.dge)')
                r('y.et <- exactTest(y.dge)')
                r('y.tp <- topTags(y.et, n = 100000)')
                r('y.pvalues <- y.et$table$PValue')
                wr = 'write.table(y.tp$table, "' + self._output + '", sep = "\t", quote = FALSE)'
                r(wr)
                finish_message = "--- edgeR with replicates is completed!"
            self._message.message_9(finish_message)
        except RRuntimeError as rre:
            self._message.message_9("Error in edgeR execution: " + str(rre))
            #raise rre