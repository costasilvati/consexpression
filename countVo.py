class CountVo(object):
    def __init__(self, inp, ann, ann_form, in_type, mod, out):
        # type: (basestring, basestring, basestring, basestring, basestring, basestring) -> object
        """
        Make a new object CountVo
        :rtype: CountVo
        :param inp:
        :param ann:
        :param ann_form:
        :param in_type:
        :param mod:
        :param out:
        count_vo = CountVo(in_count, self._exp_dao._annotation_file,
                                                  self._exp_dao._annotation_type, in_type, self._exp_dao._count_mode,
                                                  table_count)
        """
        self.command_count = "htseq-count"
        self.input_file = inp
        self.annotation_file = ann
        self.output_command = ">"
        self.output_file = out
        self.input_type = ""
        if in_type != "":
            self.input_type = "-f " + in_type
        self.annotation_form = ""
        if ann_form != "":
            self.annotation_form = "-i " + ann_form
        self.mode = ""
        if mod != "":
            self.mode = "-m " + mod

    def to_string(self):
        """

        :return:
        """

        aux = ""
        space = " "
        aux = self.command_count + space
        #aux = aux + self.input_type + space
        aux = aux + self.annotation_form + space
        aux = aux + self.mode + space
        aux = aux + self.input_file + space
        aux = aux + self.annotation_file + space
        aux = aux + self.output_command + space
        aux = aux + self.output_file + space
        print(aux)
        return aux

        # ===== TESTES DA CLASSE =====================
        # input = "~/Documents/Projeto_Juliana/Datasets/SRX080224_human_maqc-UHR1/bwa_15_12_2015-17:34/SRR299111.index.bam.sorted.bam"
        # in_type = "bam"
        # output = "output_countVo_teste.txt"
        # gtf = "~/Documents/Projeto_Juliana/Datasets/Referencias/hg19_GRCh38/hg38.gtf"
        # mode = "intersection-strict"
        # f = "gene_name"
        #
        # teste = CountVo(input, gtf, f, in_type, mode, output)
        # ok = teste.to_string()
        # print ok