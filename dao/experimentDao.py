# -*- coding: utf-8 -*-

from bo.message import Message

class ExperimentDao(object):
    """
    Object manager data of experiment
    """

    def __init__(self):
        self._message = Message()
        self._file_conf = None
        self._name_par = "NAME"
        self._replic_parm = "REPLIC"
        self._group_number_parm = "GROUP_NUMBER"
        self._group_name_parm = "GROUP_NAMES"
        self._reference_parm = "REFERENCE_GENOME"
        self._read_directory_parm = "READS_DIRECTORY"
        self._group_directory_parm = "GROUP_DIRECTORIES"
        self._paired_end_parm = "PAIRED_END"
        self._threads_parm = "THREADS"
        self._count_mode_parm = "MODE"
        self._annotation_file_parm = "ANOTATION_FILE"
        self._annotation_type_parm = "ANOTATION_TYPE"
        self._output_parm = "OUTPUT"
        self._name = ""
        self._replic = 0
        self._group_number = 0
        self._group_name = []
        self._reference = ""
        self._read_directory = ""
        self._group_directory = []
        self._paired_end = False
        self._threads = 0
        self._count_mode = ""
        self._annotation_file = ""
        self._annotation_type = ""
        self._output = ""

    def read_configuration_file(self, file):
        """
        Read file and feed class attributes, any error terminates execution
        :param file: path to config file
        :return: void
        """
        self._message.message_9("- Reading configuration file.. ----")
        conf = open(file, 'r')
        count_line = 0
        parms = {}

        for line in iter(conf):
            count_line += 1
            if line[0] != "#" and line[0] != "":
                l = line.rstrip("\n")
                p = l.split(": ")

                if p[0] in parms:
                    self._message.message_9("Parameter  " + p[0] + " is repeated!")
                else:
                    if len(p) < 2:
                        parms[p[0]] = ""
                    else:
                        parms[p[0]] = p[1]

        if self._name_par in parms:
            self._name = parms[self._name_par]
        if self._replic_parm in parms:
            self._replic = int(parms[self._replic_parm])
        if self._group_number_parm in parms:
            self._group_number = int(parms[self._group_number_parm])
        if self._group_name_parm in parms:
            self._group_name = parms[self._group_name_parm].split(',')
        if self._reference_parm in parms:
            self._reference = parms[self._reference_parm]
        if self._read_directory_parm in parms:
            self._read_directory = parms[self._read_directory_parm]
        if self._group_directory_parm in parms:
            self._group_directory = parms[self._group_directory_parm].split(',')
        if self._paired_end_parm in parms:
            self._paired_end = parms[self._paired_end_parm]
        if self._threads_parm in parms:
            self._threads = int(parms[self._threads_parm])
        if self._count_mode_parm in parms:
            self._count_mode = parms[self._count_mode_parm]
        if self._annotation_file_parm in parms:
            self._annotation_file = parms[self._annotation_file_parm]
        if self._annotation_type_parm in parms:
            self._annotation_type = parms[self._annotation_type_parm]
        if self._output_parm in parms:
            self._output = parms[self._output_parm]

# # # #================ TESTE DA CLASSE =====================================
# file = "/home/juliana/Dropbox/UTFPR/PPGBIOINFO/Projeto/RNA_tool/dao/CONFIG_tool"
# exp = ExperimentDao()
# exp.read_configuration_file(file)
# print "---"
# print exp._name