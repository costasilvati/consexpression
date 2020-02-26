from rpy2.rinterface_lib.embedded import RRuntimeError

from bo.message import Message
import rpy2.robjects as robjects
from rpy2.rinterface import *
import warnings
warnings.filterwarnings("ignore", category=RRuntimeWarning)

class Ebseq (object):

    def __init__(self, count, group, repl, out):
        """
        Inite object Ebseq
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
        self._exp_column = 1
        self._exp = "DE"

    def run_de(self, gene):
        de = 0
        if gene[self._exp_column] == self._exp:
            de = 1
        return de

    def run_ebseq(self):
        """
        Execute default analysis with EBSeq
        :return:
        """
        try:
            robjects.r('library("'+'EBSeq'+'")')
            ct = 'table <- read.csv("' + self._table_count + '",  row.names = 1, header = TRUE, stringsAsFactors=FALSE)'
            res = robjects.r(ct)
            res = robjects.r('m <- as.matrix(table)')
            grup = ""
            for ind in iter(self._groups_name):
                aux = "'" + ind + "', "
                grup = aux + grup
            grup = grup[:(len(grup) - 2)]
            # siz = 'data(m)'
            # robjects.r(siz)
            siz = 'Sizes=MedianNorm(m)'
            robjects.r(siz)
            ct = 'EBOut=EBTest(Data=m, ' \
                 'Conditions=as.factor(rep(' \
                'c(' + grup + '),each=' + str(self._replic) + ')), sizeFactors=Sizes, maxround=5)'
            robjects.r(ct)
            ct = 'EBDERes=GetDEResults(EBOut, FDR=0.05)'
            robjects.r(ct)
            wr = 'write.table(EBDERes$Status, file="' + self._output + '", sep = "\t", quote = FALSE)'
            robjects.r(wr)
            self._message.message_9("--- EBSeq: is completed!")
        except RRuntimeError as rre:
            self._message.message_9("Error in baySeq execution: " + str(rre))
            raise rre
