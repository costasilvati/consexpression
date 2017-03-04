# -*- coding: utf-8 -*-

from bo.samseq import SamSeq
from bo.limmavoom import LimmaVoom
from bo.baySeq import BaySeq
from bo.deseq import DESeq
from bo.ebseq import Ebseq
from bo.edgeR import EdgeR
from bo.noiseq import Noiseq

class Consensus (object):
    """
    Make a consensus of expression analysis results
    """

    def __init__(self, expression_results):

        self._exp_results = expression_results
        self._likelihood = 0.95
        self._fdr = 0.1
        self._lfc = 2
        self._pvalue = 0.05
        self._qvalue = 1

    def run_consensus(self):
        """
        Make a consensus of expression analysis results
        :return:
        """
        for method in self._exp_results:
