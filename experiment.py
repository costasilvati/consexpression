# -*- coding: utf-8 -*-
import os
import sys
import glob
from bo.message import Message
from dao.experimentDao import ExperimentDao
from vo.mappVo import MappVo
from bo.mappBo import MappBo
from vo.countVo import CountVo
from bo.countBo import CountBo
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
    _count: CountBo

    def __init__(self):

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



    def init_experiment(self, exp, file):
        """
        Iniatialize experiment
        :param exp:
        :param file: config file
        :return:
        """
        assert isinstance(exp, ExperimentDao)
        self._exp_dao = exp
        self._exp_dao.read_configuration_file(file)
        self._exp_dao._name = self.name_valid(self._exp_dao._name)
        self.rep_valid(self._exp_dao._replic)
        self.group_number_valid(self._exp_dao._group_number)
        ref = self._exp_dao._reference

        if self._exp_dao._reference == "":
            print ("You don't have a refserence genome... Expression analyse need a table count with mapping reads")
            self._merged_table_out = input("Type absolute path to table count")
        elif ref != "" and (self.extension_valid(ref, "fa") or self.extension_valid(ref, "fasta")):
            self._reference = self._exp_dao._reference # == ref
        else:
             self.message.message_3("REFERENCE FILE ")
             exit(0)

        self.directory_valid(self._exp_dao._read_directory, "reads")

        for i in iter(self._exp_dao._group_directory):
            reads = self._exp_dao._read_directory + "/" + str(i)
            self.directory_valid(reads, "group")
            # Get fastq reads
            path_find = self._exp_dao._read_directory + "/" + i + "/"
            self._fastq.append(self.get_reads_file(path_find))


        if self._exp_dao._paired_end == True:
            self.message.message_8("The sequence is paired-end. CONSEXPRESSION dont make paired-end analysis")
            exit(0)
        else:
            self.message.message_8("The sequence is single-end")

    def name_valid(self, name):
        """
        Verify the name: if is empty change to default name
        :param name: name of experiment
        :return: boolean
        """
        if len(name) == 0:
            name = "consexpression"
            self.message.message_7("Experiment name is empty! The name was changed to consexpression")
        return name

    def rep_valid(self, rep):
        """
        Verify if number of replicates technical and biological is valid (>= 1).
        basestring :param rep_t:
        basestring :param rep_b:
        void :return:
        """
        ok = False
        if rep >= 1:
            self.message.message_1("replics")
        else:
            self.message.message_2("1 replic or more (technique or biological)")
            self.message.message_3("number of replics in line 5 - 6")
            exit()

    def extension_valid(self, path, extension):
        """
        Verify if the extension file is the expected
        :param path:
        :param extension:
        :return: boolean
        """
        var_ret = False
        if str.endswith(path, extension):
            var_ret = True
        else:
            var_ret = False

        return var_ret

    def directory_valid(self, path, type):
        """
        Verify if path is a directory
        :param path: path of file
        :param type: file is reference genome, reads?
        :return: void
        """
        ok = os.path.isdir(path)
        if ok:
            self.message.message_1("directory " + type + ": " + path)
        else:
            self.message.message_2("a valid directory " + type + " path")
            self.message.message_3(" the directory " + type + " path (line 9)")
            exit()

    def group_number_valid(self, group_n):
        """
        Verify the number of groups. The minimal is one
        int :param group_n:
        void :return:
        """
        assert isinstance(group_n, int)

        if group_n >= 1:
           self.message.message_1("group number")
        else:
            self.message.message_2("1 group or more.")
            self.message.message_3("the number of gruoups in line 7")
            exit()

    def file_valid(self, path):
        """
        Verify if path is a file
        basestring :param path:
        void :return:
        """
        if os.path.isfile(path):
            self.message.message_1(" file: " + path)
        else:
            self.message.message_2(" a valid reference file")
            self.message.message_3(" the reference file path (line 7)")
            exit()

    def exceute_mapp_count(self):
        """
        Execute Tophat and htseq-count
        :return:
        """
        ref = self._reference
        thread = self._exp_dao._threads
        sing = self._exp_dao._paired_end
        n = self._exp_dao._name
        path_find = []
        for grp in iter(self._fastq):
            for grp_file in iter(grp):
                bar = 1 + grp_file.rfind('/')
                out_mapp = grp_file[:bar] + n + "/" + grp_file[bar:]
                dir = grp_file[:bar] + n
                if os.path.isdir(dir):
                    pass
                else:
                    os.mkdir(dir, 0o755)
                out_mapp = out_mapp.replace('fastq', 'sam')
                path_find.append(out_mapp)
                mapp_vo = MappVo("TopHat", ref, grp_file, "", thread, out_mapp, "", sing)
                self._mapp_bo = MappBo(mapp_vo)
                mapp_exe = self._mapp_bo.execute_mapp()
                if mapp_exe == 0:
                    dot = out_mapp.rfind('.')
                    in_type = out_mapp[dot + 1:]
                    bar = 1 + out_mapp.rfind('/')
                    table_count = out_mapp[bar:dot]
                    table_count = out_mapp[:bar] + table_count + "_table_count.txt"
                    self._count_table.append(table_count)
                    in_count = out_mapp + "/accepted_hits.sam"

                    count_vo = CountVo(in_count, self._exp_dao._annotation_file,
                                                  self._exp_dao._annotation_type, in_type, self._exp_dao._count_mode,
                                                  table_count)
                    self._count = CountBo(count_vo)
                    if self._count.execute_count() == 0:
                        self.message.message_8("Count Sucsessfull!!!")
                    else:
                        self.message.message_4("Error in counting mapped reads...")
                else:
                    self.message.message_4("Task: Mapping don't run correctly.")
            self._out_mapp.append(path_find)


    def get_reads_file(self, dir):
        """
        Get all fastq path of dir
        :param dir: path to folder of fastq sample
        :return: array of reads file path
        """
        fastq_file = []
        path = dir + "*.fastq" #serach
        for file in glob.glob(path):
            fastq_file.append(file)
        if len(fastq_file) == 0:
            self.message.message_7("*Not found files FASTQ in directorie " + dir)
        return fastq_file

    def execute_expression_analysis(self):
        """
        Make analysis with counts data for mapping
        :return:
        """
        print ("Expression analisys start...")
        n = "consexpression"
        out_merge_table = ""
        if self._exp_dao._reference != "":
            out_merge_table = self._exp_dao._read_directory + "/" + self._exp_dao._name + "_table_count.txt"
            self.execute_merge_table(self._count_table, out_merge_table)
        else:
            out_merge_table = self._merged_table_out
        # 1 ------------------ edgeR -----------------
        print("---- edgeR START! ------------")
        out_expression = self._exp_dao._output +"/"+ self._exp_dao._name
        out_edger =  out_expression +"_edger.csv"
        self._edger = EdgeR(out_merge_table, self._exp_dao._group_name, self._exp_dao._replic, out_edger)
        self._edger.run_edger()
        # 2 ------------- BaySeq --------------------
        print("---- baySeq START! ------------")
        out_bayseq =  out_expression + "_baySeq.csv"
        self._bayseq = BaySeq(out_merge_table, self._exp_dao._group_name, self._exp_dao._replic, out_bayseq)
        self._bayseq.run_bayseq()
        # 3 ------------- DESeq --------------------
        print("---- DESeq START! ------------")
        out_deseq =  out_expression + "_DESeq.csv"
        self._deseq = DESeq(out_merge_table, self._exp_dao._group_name, self._exp_dao._replic, out_deseq)
        self._deseq.run_deseq()
        # 4 ------------- NOISeq --------------------
        print("---- NOISeq START! ------------")
        out_noiseq =  out_expression + "_NOISeq.csv"
        self._noiseq = Noiseq(out_merge_table, self._exp_dao._group_name, self._exp_dao._replic, out_noiseq)
        self._noiseq.run_noiseq()
        # 5 ------------- EBSeq --------------------
        print("---- EBSeq START! ------------")
        out_ebseq =  out_expression + "_EBSeq.csv"
        self._ebseq = Ebseq(out_merge_table, self._exp_dao._group_name, self._exp_dao._replic, out_ebseq)
        self._ebseq.run_ebseq()
        # 6 ------------- SAMSeq --------------------
        print("---- SAMSeq START! ------------")
        # out_samseq =  out_expression + "_SAMSeq.csv"
        # self._samseq = SamSeq(out_merge_table, self._exp_dao._group_name, self._exp_dao._replic, out_samseq)
        # self._samseq.run_samseq()
        # 7 ------------- limma-voom --------------------
        print("---- limma START! ------------")
        out_limmavoom =  out_expression + "_limmavoom.csv"
        self._limmavoom = LimmaVoom(out_merge_table, self._exp_dao._group_name, self._exp_dao._replic, out_limmavoom)
        self._limmavoom.run_limmavoom()

    def execute_conseus(self, out):
        gene_de = {}
        read_bay = open(self._bayseq._output, 'r')
        c_b = 0
        for line in iter(read_bay):
            if c_b > 0:
                gene = line.split("\t")
                v = self._bayseq.run_de(gene)
                if gene[0] in gene_de:
                    aux = gene_de[gene[0]]
                    gene_de[gene[0]] = int(aux) + int(v)
                else:
                    gene_de[gene[0]] = int(v)
            c_b += 1
        read_bay.close()

        # ---- edger
        read_edger = open(self._edger._output, 'r')
        c_b = 0
        for line in iter(read_edger):
            if c_b > 0:
                gene = line.split("\t")
                v = self._edger.run_de(gene)
                if gene[0] in gene_de:
                    aux = gene_de[gene[0]]
                    gene_de[gene[0]] = int(aux) + int(v)
                else:
                    gene_de[gene[0]] = int(v)
            c_b += 1
        read_edger.close()

        #--- deseq
        read_deseq = open(self._deseq._output, 'r')
        c_b = 0
        for line in iter(read_deseq):
            if c_b > 0:
                gene = line.split("\t")
                v = self._deseq.run_de(gene)
                if gene[1] in gene_de:
                    aux = gene_de[gene[1]]
                    gene_de[gene[1]] = int(aux) + int(v)
                else:
                    gene_de[gene[1]] = int(v)
            c_b += 1
        read_deseq.close()

        # --- noiseq
        read_noiseq = open(self._noiseq._output, 'r')
        c_b = 0
        for line in iter(read_noiseq):
            if c_b > 0:
                gene = line.split(",")
                v = self._noiseq.run_de(gene)
                if gene[0] in gene_de:
                    aux = gene_de[gene[0]]
                    gene_de[gene[0]] = int(aux) + int(v)
                else:
                    gene_de[gene[0]] = int(v)
            c_b += 1
        read_noiseq.close()

        # --- samseq
        if self._samseq is None:
            print("SAMSeq results not found")
        else:
            read_samseq = open(self._samseq._output, 'r')
            c_b = 0
            for line in iter(read_samseq):
                if c_b > 0:
                    gene = line.split("\t")
                    v = self._samseq.run_de(gene)
                    if gene[1] in gene_de:
                        aux = gene_de[gene[1]]
                        gene_de[gene[1]] = int(aux) + int(v)
                    else:
                        gene_de[gene[1]] = int(v)
                c_b += 1
            read_samseq.close()

        # --- limma
        if self._exp_dao._replic >= 2:
            read_limma = open(self._limmavoom._output, 'r')
            c_b = 0
            for line in iter(read_limma):
                if c_b > 0:
                    gene = line.split("\t")
                    v = self._limmavoom.run_de(gene)
                    if gene[0] in gene_de:
                        aux = gene_de[gene[0]]
                        gene_de[gene[0]] = int(aux) + int(v)
                    else:
                        gene_de[gene[0]] = int(v)
                c_b += 1
            read_limma.close()
        else:
            print("limma require more than one replics")

        # --- ebseq
        read_ebseq = open(self._ebseq._output, 'r')
        c_b = 0
        for line in iter(read_ebseq):
            if c_b > 0:
                gene = line.split("\t")
                v = self._ebseq.run_de(gene)
                if gene[0] in gene_de:
                    aux = gene_de[gene[0]]
                    gene_de[gene[0]] = int(aux) + int(v)
                else:
                    gene_de[gene[0]] = int(v)
            c_b += 1
        read_ebseq.close()

        #--- write results
        header = 'gene, indications'
        out_cons = open(out, 'w')
        out_cons.write(header)
        names = gene_de.keys()
        print(len(names))
        for i in iter(names):
            if (gene_de[i]) >= 4:
                out_cons.write("\n" + i + "," + str(gene_de[i]))

    def execute_merge_table(self, out_mapp_list, out_name):
        """
        Make a merge table whit counts
        :param out_mapp_list:
        :param out_name:
        :return:
        """
        self._count.merge_table_count(out_mapp_list, out_name, self._exp_dao._group_name)



if __name__ == "__main__":
    expBo = Experiment()
    config_file = "/Users/julianacostasilva/PycharmProjects/consexpression/dao/CONFIG_tool.txt" # sys.argv[1]
    exp_d = ExperimentDao()
    expBo.init_experiment(exp_d, config_file)
    if(exp_d._reference != ""):
        expBo.exceute_mapp_count()
    expBo.execute_expression_analysis()
    expBo.execute_conseus(exp_d._output+'/consensus.txt')