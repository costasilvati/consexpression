# -*- coding: utf-8 -*-
import os
import sys
import glob
pmName = input('bo.message')
pm = __import__(pmName)
from bo.message import Message
from dao.experimentDao import ExperimentDao
from bo.edgeR import EdgeR
from bo.baySeq import BaySeq
from bo.deseq import DESeq
from bo.noiseq import Noiseq
from bo.ebseq import Ebseq
from bo.limmavoom import LimmaVoom
from bo.samseq import SamSeq

class Experiment(object):
    """
        Business object of Experiment
    """
    def __init__(self):
        print("-----------------------------------"+
            "\n -- Welcome to consExpression -- "+
            "\n---------------------------------\n")
        self._exp_dao = None
        self._reference = None
        self._transcript = False
        self._count = None
        self._expression = None
        self._mapp_bo = None
        self.message = Message()
        self._fastq = []
        self._out_mapp = []
        self._count_table = []
        self._merged_table_out = None
        self._edger = None
        self._bayseq = None
        self._deseq = None
        self._noiseq = None
        self._ebseq = None
        self._samseq = None
        self._limmavoom = None


    def execute_expression_analysis(self):
        """
        Make analysis with counts data for mapping
        :return:
        """
        print ("Expression analisys start...")
        n = "consexpression"
        out_merge_table = ''
        self.execute_merge_table(self._count_table, out_merge_table)
        # 1 ------------------ edgeR -----------------
        out_edger = self._exp_dao._read_directory + "/" + self._exp_dao._name + "_edger.csv"
        self._edger = EdgeR(out_merge_table, self._exp_dao._group_name, self._exp_dao._replic, out_edger)
        self._edger.run_edger()
        # 2 ------------- BaySeq --------------------
        out_bayseq = self._exp_dao._read_directory + "/" + self._exp_dao._name + "_baySeq.csv"
        self._bayseq = BaySeq(out_merge_table, self._exp_dao._group_name, self._exp_dao._replic, out_bayseq)
        self._bayseq.run_bayseq()
        # 3 ------------- DESeq --------------------
        out_deseq = self._exp_dao._read_directory + "/" + self._exp_dao._name + "_DESeq.csv"
        self._deseq = DESeq(out_merge_table, self._exp_dao._group_name, self._exp_dao._replic, out_deseq)
        self._deseq.run_deseq()
        # 4 ------------- NOISeq --------------------
        out_noiseq = self._exp_dao._read_directory + "/" + self._exp_dao._name + "_NOISeq.csv"
        self._noiseq = Noiseq(out_merge_table, self._exp_dao._group_name, self._exp_dao._replic, out_noiseq)
        self._noiseq.run_noiseq()
        # 5 ------------- EBSeq --------------------
        out_ebseq = self._exp_dao._read_directory + "/" + self._exp_dao._name + "_EBSeq.csv"
        self._ebseq = Ebseq(out_merge_table, self._exp_dao._group_name, self._exp_dao._replic, out_ebseq)
        self._ebseq.run_ebseq()
        # 6 ------------- SAMSeq --------------------
        out_samseq = self._exp_dao._read_directory + "/" + self._exp_dao._name + "_SAMSeq.csv"
        self._samseq = SamSeq(out_merge_table, self._exp_dao._group_name, self._exp_dao._replic, out_samseq)
        self._samseq.run_samseq()
        # 7 ------------- limma-voom --------------------
        out_limmavoom = self._exp_dao._read_directory + "/" + self._exp_dao._name + "_limmavoom.csv"
        self._limmavoom = LimmaVoom(out_merge_table, self._exp_dao._group_name, self._exp_dao._replic, out_limmavoom)
        self._limmavoom.run_limmavoom()

if __name__ == "__main__":
    expBo = Experiment()
    config_file = sys.argv[1]
    exp_d = ExperimentDao()
    expBo.init_experiment(exp_d, config_file)
    expBo.exceute_mapp_count()
    expBo.execute_expression_analysis()
#    expBo.execute_conseus('consensus_teste.txt')
