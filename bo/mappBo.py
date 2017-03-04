# coding=utf-8

from vo.mappVo import MappVo
from bo.message import Message
import multiprocessing
import glob
import subprocess
import os


class MappBo(object):
    """
    This class make rules of validate information and command, to execute Mapp tools
    """

    def __init__(self, mapp):
        assert isinstance(mapp, MappVo)
        self._map_vo = mapp
        self._reads_file = []
        self.message = Message()

    def threads_conf(self,threads_vo):
        """
        Alter threads larger to default of system
        :param threads_vo: number of threads
        :return: void
        """
        threads_sys = multiprocessing.cpu_count()
        if threads_vo < threads_sys:
            self._map_vo._threads_value = threads_sys-1
            self.message.message_9("The threads nunber defined is " +
                                   str(threads_vo) +
                                   ", but the system have only "
                                   + str(threads_sys))
            self.message.message_9("---> Number of threads was change to "
                                   + str(threads_sys-1))
        self.message.message_9("Successful! Threads configuration is ok!")

    def execute_mapp(self):
        """
        Execute the command: 0 is ok, 1 is fail mapped task
        :return: int
        """
        self.threads_conf(self._map_vo._threads_value)
        n = self.make_bowtie2_index(self._map_vo._index_name)
        self._map_vo._index_name = n
        text = self._map_vo.to_string()
        return_code = subprocess.call(text, shell=True)
        return return_code

    def make_bowtie2_index(self, index):
        """
        Execute command to make a bowtie2 index if do not exists
        :param index: fasta file reference to mapp
        :return: name of generated index
        """
        dot = index.rfind('.f')
        name = index[:dot]
        if os.path.isfile(name + ".1.bt2"):
            return name
        else:
            command = "bowtie2-build " + index + " " + name
            if subprocess.call(command, shell=True) == 0:
                return name
            else:
                self.message.message_4("Error in index build")
                return ""

# #===== TESTES DA CLASSE =====================
# name = "Bowtie2"
# index_name = "/home/juliana/Documents/Projeto_Juliana/Datasets/Referencias/GRCh38.p5/GCA_000001405.20_GRCh38.p5_genomic.fna"
# threads_value = 3
# reads1_name = "/home/juliana/Documents/Eliandro-UEL/E1_S1_L001_R1_001_prinseq_1.fastq"
# reads2_name = "/home/juliana/Documents/Eliandro-UEL/E1_S1_L001_R2_001_prinseq_2.fastq"
# output_name = "/home/juliana/Documents/Testes_RNATool/eliandro_uel.sam"
# map = vo.MappVo.MappVo(name,index_name,reads1_name, reads2_name, threads_value,output_name,"",False)
# mapbo = MappBo(map)
# mapbo.make_bowtie2_index(index_name)
# # map.parm_mapp()
# # teste = map.to_string()
# # print teste
# # print teste