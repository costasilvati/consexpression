
from bo.message import Message
import rpy2.robjects as robjects
from rpy2.rinterface import *
import warnings
warnings.filterwarnings("ignore", category=RRuntimeWarning)

class SamSeq (object):

    def __init__(self, count, group, repl, out):
        """
        Inite object Ebseq
        :param count:
        :param group:
        :param repl:
        :param out:
        """
        robjects.r['options'](warn=-1)
        self._table_count = count
        self._groups_name = group
        self._replic = repl
        self._output = out
        self._class = '"Two class unpaired"'
        self._message = Message()
        self._fd_column = 4
        self._qvalue_column = 5
        self._qvalue = 1
        self._fd = 2

    def run_de(self, gene):
        de = 0
        fd = float(gene[self._fd_column])
        qv = float(gene[self._qvalue_column])
        if fd <= self._fd and fd <= self._qvalue:
            de = 1
        return de


    def run_samseq(self):
        """
        Execute default analysis with SAMSeq
        :return:
        """
        try:
            if len(self._groups_name) > 2:
                self._class = '"Multiclass"'

            robjects.r('library("'+'samr'+'")')
            res = robjects.r('table <- read.csv("' + self._table_count + '",  row.names = 1, header = TRUE, stringsAsFactors=FALSE, sep = "' + ',' + '")')
            res = robjects.r('m <- as.matrix(table)')

            grup = ""
            for ind in iter(self._groups_name):
                grup = grup + '"' + ind + '",'
            grup = grup[:(len(grup) - 1)]

            cm = 'SAMseq.test = SAMseq(m, as.factor(rep(c('
            cm = cm + grup + '),each=' + str(self._replic) + ')), resp.type = '+ self._class + ', geneid = rownames(m), genenames = rownames(m), nperms = 100)'
            #print(cm)
            res = robjects.r(cm)
            res = robjects.r('SAMseq.result.table = rbind(SAMseq.test$siggenes.table$genes.up, SAMseq.test$siggenes.table$genes.lo)')
            res = robjects.r('SAMseq.score = rep(0, nrow(m))')
            res = robjects.r('SAMseq.score[match(SAMseq.result.table[,1], rownames(m))] = as.numeric(SAMseq.result.table[,3])')
            res = robjects.r('SAMseq.FDR = rep(1, nrow(m))')
            res = robjects.r('SAMseq.FDR[match(SAMseq.result.table[,1], rownames(m))] = as.numeric(SAMseq.result.table[,5])/100')
            wr = 'write.table(SAMseq.result.table, file="' + self._output + '", sep = "\t", quote = FALSE)'
            robjects.r(wr)
            self._message.message_9("--- SAMSeq: is completed!")
        except RRuntimeError as rre:
            self._message.message_9("Error in SAMSeq execution: " + str(rre))
            # raise rre