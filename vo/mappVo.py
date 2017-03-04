from bo.message import Message


class MappVo(object):
    """
    Record values to run Mapp methos
    """

    #mapp_vo = MappVo("TopHat", ref, grp_file, "", thread, out_mapp, "", sing)
    def __init__(self, name, index, read1_n, read2_n, threads, out, other, single_end):
        """
        Meka a first construction of object with param by user
        :param name:
        :param index:
        :param reads_dir:
        :param threads:
        :param out:
        :param other:
        :param single_end:
        """
        self._index_parm = ""
        self._reads1_parm = ""
        self._reads2_parm = ""
        self._threads_parm = ""
        self._output_parm = ""
        self._command_parm = ""
        self._sep = " "
        self._name = name
        #nIndex = index[(1 + index.rfind('/')):] # (index.rfind('.'))
        self._index_name = index # nIndex
        self._reads1_name = read1_n
        self._reads2_name = read2_n
        self._threads_value = threads
        self._output_name = out
        self._output_type = "--no-convert-bam "
        self._other_conf = other
        self._paired_end = single_end
        self._message = Message()
        self.parm_mapp()

    def parm_mapp(self):
        """
        Make command parameter by mapping tool
        """

        if self._name == "BWA":
            self._command_parm = "bwa mem "
            self._threads_parm = "-t "
            self._output_parm = "> "

        elif self._name == "Bowtie2":
            self._command_parm = "bowtie2 "
            self._index_parm = "-x "
            self._threads_parm = "-p "
            self._output_parm = "-S "

            if self._paired_end == 'True':
                self._reads1_parm = "-1 "
                self._reads2_parm = "-2 "
            else:
                self._reads1_parm = "-U "

        elif self._name == "TopHat":
            self._command_parm = "tophat2 "
            self._threads_parm = "-p "
            self._output_parm = "--output-dir "
        else:
            self._message.message_4("Mapping " + self._name + " not found!")
            exit()

    def to_string(self):
        """
        Return a command, this command used to run a Mapping tool
        str :return:
        """
        aux = ""
        if self._name == "Bowtie2":
            aux = self._command_parm + self._index_parm + self._index_name + self._sep
            aux = aux + self._threads_parm + str(self._threads_value) + self._sep
            aux = aux + self._reads1_parm + self._reads1_name + self._sep
            if self._paired_end == 'True':
                aux = aux + self._reads2_parm + self._reads2_name + self._sep
            aux = aux + self._output_parm + self._output_name + self._sep
            aux = aux + self._other_conf
            return aux
        elif self._name == "BWA":
            aux = self._command_parm + self._index_parm + self._index_name + self._sep
            aux = aux + self._threads_parm + str(self._threads_value) + self._sep
            aux = aux + self._reads1_parm + self._reads1_name + self._sep
            if self._paired_end == 'True':
                aux = aux + self._reads2_parm + self._reads2_name + self._sep
            aux = aux + self._output_parm + self._output_name + self._sep
            aux = aux + self._other_conf
            return aux
        elif self._name == "TopHat":
            aux = self._command_parm
            aux = aux + self._threads_parm + str(self._threads_value) + self._sep
            aux = aux + self._output_type
            aux = aux + self._other_conf
            aux = aux + self._output_parm + self._output_name + self._sep
            #print(aux + self._index_name)
            aux = aux + self._index_parm + self._index_name + self._sep
            aux = aux + self._reads1_parm + self._reads1_name + self._sep
            if self._paired_end == 'True':
                aux = aux + self._reads2_parm + self._reads2_name + self._sep

            return aux
        else:
            return aux