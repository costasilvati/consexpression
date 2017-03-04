from vo.countVo import CountVo
from bo.message import Message
import subprocess

class CountBo(object):
    """
    This object define business rules to make count table execution
    """

    def __init__(self, count):
        """
        Test the doc of constructor class
        :param count:
        """
        assert isinstance(count, CountVo)
        self._counter = count
        self.message = Message()


    def annotation_format(self):
        """
        Verify format of annotation file (default: GTF | GFF)
        :return: void
        """
        bar = self._counter.annotation_file.rfind('.')
        name = self._counter.annotation_file[bar:]

        if name != 'gtf' and name != 'gff':
            self.message.message_4('File extension of annotation file can be only GTF or GFF.')

    def execute_count(self):
        """
        Execute command htseq-count
        :return: int subprocess
        """
        text = self._counter.to_string()
        return_code = subprocess.call(text, shell=True)
        return return_code

    def merge_table_count(self, list_file, out, groups_name):
        """
        Make a table whit count of all samples
        :param list_file: array count files
        :param out: text file line (gene) column (sample) data (count mapped)
        :param groups_name: treatment of samples
        :return:
        """
        n_g = len(groups_name)
        group_count = 0
        rep = int(len(list_file) / n_g)
        rep_count = 1
        out_file = None
        gene = {}
        no_genes = {'__no_feature': 0, '__ambiguous': 0, '__too_low_aQual': 0, '__not_aligned': 0,'__alignment_not_unique': 0,
                    'not_aligned':0, 'no_feature':0, 'ambiguous':0, 'too_low_aQual':0, 'alignment_not_unique':0}
        out_file = open(out, 'w')
        out_file.write("gene")
        # loop table count by samples
        for ind in iter(list_file):
            op = open(ind, 'r')
            if rep_count <= rep:
                out_file.write("," + groups_name[group_count] + str(rep_count))
            else:
                rep_count = 1
                group_count += 1
                out_file.write("," + groups_name[group_count]+str(rep_count))

            for line in iter(op):
                line = line.rstrip()
                text = line.split("\t")
                if text[0] in no_genes:
                    pass
                else:
                    if text[0] in gene:
                        aux = gene[text[0]]
                        aux = aux + ',' + text[1]
                        gene[text[0]] = aux
                    else:
                        gene[text[0]] = text[1]
            op.close()
            rep_count += 1
        names = gene.keys()
        for i in iter(names):
            out_file.write("\n" + i + "," + str(gene[i]))
        out_file.close()