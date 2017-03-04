
from bo.message import Message
import rpy2.robjects as robjects
from rpy2.rinterface import *
import warnings
warnings.filterwarnings("ignore", category=RRuntimeWarning)

class LimmaVoom (object):

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
        self._logfc_column = 2
        self._pvalue_column = 5
        self._logfc = 2
        self._pvalue = 0.05

    def run_de(self, gene):
        de = 0
        lfc = float(gene[self._logfc_column])
        pv = float(gene[self._pvalue_column])
        if lfc >= self._logfc or lfc <= -self._logfc:
            if pv >= self._pvalue:
                de = 1
        return de

    def run_limmavoom(self):
        """
        Execute default analysis with Limma-voom
        :return:
        """
        if self._replic == 1:
            self._message.message_4("limma-voom require more than one replics.")
            self._message.message_9("--- limma-voom: is kipped!")
        else:
            try:
                robjects.r('library("'+'edgeR'+'")')
                robjects.r('library("' + 'limma' + '")')
                ct = 'table <- read.csv("' + self._table_count + '",  row.names = 1, header = TRUE, stringsAsFactors=FALSE)'
                res = robjects.r(ct)
                res = robjects.r('m <- as.matrix(table)')
                res = robjects.r('nf = calcNormFactors(m, method = "TMM")')

                grup = ""
                for ind in iter(self._groups_name):
                    grup = grup + ('"' + ind + '",') * self._replic
                grup = grup[:(len(grup) - 1)]
                robjects.r('condition = factor(c(' + grup + '))')

                res = robjects.r('voom.data <- voom(m, model.matrix(~factor(condition)), lib.ize = colSums(m) * nf)')
                res = robjects.r('voom.data$genes = rownames(m)')
                res = robjects.r('voom.fitlimma = lmFit(voom.data, design=model.matrix(~factor(condition)))')
                res = robjects.r('voom.fitbayes = eBayes(voom.fitlimma)')
                res = robjects.r('voom.pvalues = voom.fitbayes$p.value[, 2]')
                res = robjects.r('voom.adjpvalues = p.adjust(voom.pvalues, method="BH")')
                var = 'design <- c(' + '1,' * self._replic + '2,'*self._replic
                var = var[:(len(var) - 1)] + ')'
                res = robjects.r(var)
                res = robjects.r('data <- topTable(voom.fitbayes, coef=ncol(design), number=1000000)')
                wr = 'write.table(data, file="' + self._output + '", sep = "\t", quote = FALSE)'
                robjects.r(wr)
                self._message.message_9("--- limma-voom: is completed!")
            except RRuntimeError as rre:
                self._message.message_9("Error in limma-voom execution: " + str(rre))
                # raise rre