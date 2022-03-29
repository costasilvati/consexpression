# coding=utf-8
from rpy2.rinterface_lib.embedded import RRuntimeError

from bo.message import Message
import rpy2.robjects as robjects
from rpy2.rinterface import *
import warnings
warnings.filterwarnings("ignore", category=RRuntimeWarning)


class BaySeq(object):
    """

    Commands to run BaySeq expression analysis
    """

    def __init__(self, count, group, repl, output):
        """
        Define the edgeR object
        :param count:
        :param group:
        :param repl:
        :param output:
        """
        self._table_count = count
        self._groups_name = group
        self._replic = repl
        self._output = output
        self._message = Message()
        self._likelihood_column = 2 + len(group)*repl
        self._fdr_de_column = 4 + len(group)*repl
        self._likelihood = 0.95
        self._fdr = 0.1

    def run_de(self, gene):
        de = 0
        try:
            fdr = float(gene[self._fdr_de_column])
            like = float(gene[self._likelihood_column])
            if fdr <= self._fdr and like > self._likelihood:
                de = 1
        except ValueError:
            de = 0
        return de

    def run_bayseq(self):
        """
        Execute default analysis with baySeq
        :return:
        """
        try:
            res = robjects.r('library("parallel")')
            res = robjects.r('library("stats4")')
            res = robjects.r('library("BiocGenerics")')
            res = robjects.r('library("S4Vectors")')
            res = robjects.r('library("IRanges")')
            res = robjects.r('library("GenomeInfoDb")')
            res = robjects.r('library("abind")')
            # res = robjects.r('library("perm")')
            res = robjects.r('library("GenomicRanges")')
            res = robjects.r('library("baySeq")')

            res = robjects.r('if(require("parallel")) cl <- makeCluster(4) else cl <- NUL')
            ct = 'table <- read.csv("' + self._table_count + '",  row.names = 1, header = TRUE, stringsAsFactors = FALSE)'
            res = robjects.r(ct)
            res = robjects.r('m <- as.matrix(table)')
            replicates = ""
            assert isinstance(self._replic, int)
            for ind in iter(self._groups_name):
                aux = "'" + ind + "', "
                replicates = replicates + aux * self._replic
            replicates = replicates[:(len(replicates) - 2)]
            replicates = 'replicates <- c(' + replicates + ')'
            res = robjects.r(replicates)
            groups = 'groups <- list(NDE = c('+ "1," * len(self._groups_name)
            groups = groups[:(len(groups) - 1)] + ')'
            groups = groups + ', DE = c('+ '1,' * self._replic
            groups = groups + '2,' * self._replic
            groups = groups[:(len(groups) - 1)] + "))"
            res = robjects.r(groups)
            res = robjects.r('CD <- new("countData", data = m, replicates = replicates, groups = groups)')
            res = robjects.r('libsizes(CD) <- getLibsizes(CD)')
            res = robjects.r('CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl, equalDispersions = TRUE)')
            res = robjects.r('CD <- getLikelihoods(CD, prs=c(0.5, 0.5), pET="BIC", cl=cl)')
            # CD.posteriors.DE < - exp(CD @ posteriors)[, 2]
            res = robjects.r('write.table(topCounts(CD, group = "DE", number = 65000, normaliseData = TRUE), "' + self._output +'", sep="\t", quote = FALSE)')
            self._message.message_9("--- baySeq is completed!")
        except RRuntimeError as rre:
            self._message.message_9("Error in baySeq execution: " + str(rre))
            raise rre

#========================= TESTE da CLASSE==============
# inp = '/Volumes/SD128/bioconvergencia/reads_RNApa/kallisto_quant_scRNA_apa_tpm_tab.csv'
# gr = ["0b", "pb"]
# rp = 2
# out = 'RNApa_apa_1B_0B-consexpression_deseq.csv'
# b = BaySeq(inp, gr, rp, out)
# b.run_bayseq()
# Error in baySeq execution: Error in file(file, "rt") : não é possível abrir a conexão
# read_bay = open(inp, 'r')
# c_b = 1
# for line in iter(read_bay):
# #     #print('--' + line)
# #     if c_b > 0:
# #        gene = line.split("\t")
# #        print(gene[0])
# #        v = b.run_de(gene)
# #        print('--> '+ str(v))
# #     c_b += 1