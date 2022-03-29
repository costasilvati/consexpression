htseq-count -f sam -r pos -t exon -i gene_id --additional-attr gene_name -c L2_ScXa_i_10_ntd_count.csv -n 6 L2_ScXa_i_10_ntd.sam ../../marcelozerillo/SP80-3280_MMZ.gff3

htseq-count -f sam -r pos -t exon --additional-attr gene_name -c L2_ScXa_i_10_ntd_count.csv -n 6 L2_ScXa_i_10_ntd.sam ../../marcelozerillo/SP80-3280_MMZ.gff3
htseq-count -f sam -r pos -t gene --additional-attr gene_name -c L2_ScXa_i_10_ntd_count.csv -n 6 L2_ScXa_i_10_ntd.sam ../../marcelozerillo/SP80-3280_MMZ.gff3
htseq-count -f sam -r pos -t gene -i ID --additional-attr gene_name -c L2_ScXa_i_10_ntd_count.csv -n 6 L2_ScXa_i_10_ntd.sam ../../marcelozerillo/SP80-3280_MMZ.gff3
#feature evm.model.QPEU01000002.1.1.exon1 does not contain a 'Name' attribute
htseq-count -f sam -r pos -t exon -i Name --additional-attr gene_name --additional_attr exon_number -c L2_ScXa_i_10_ntd_count.csv -n 6 L2_ScXa_i_10_ntd.sam ../../marcelozerillo/SP80-3280_MMZ.gff3
htseq-count -f sam -r pos -t exon -i Name --additional-attr gene_name --additional_attr exon_number -c L2_ScXa_i_10_ntd_count.csv -n 6 L2_ScXa_i_10_ntd.sam ../../marcelozerillo/SP80-3280_MMZ.gff3
# Feature evm.model.QPEU01000002.1.1.exon1 does not contain a 'Name' attribute

htseq-count -f sam -r pos -t exon -i name --additional-attr gene_name --additional_attr exon_number -c L2_ScXa_i_10_ntd_count.csv -n 6 L2_ScXa_i_10_ntd.sam ../../marcelozerillo/SP80-3280_MMZ.gff3

htseq-count -f sam -r pos -t exon -i name --additional-attr gene_name -c L2_ScXa_i_10_ntd_count.csv -n 6 L2_ScXa_i_10_ntd.sam ../../marcelozerillo/SP80-3280_MMZ.gff3

htseq-count -f sam -r pos -t exon -i sequence-region --additional-attr Name -c L2_ScXa_i_10_ntd_count.csv -n 6 L2_ScXa_i_10_ntd.sam ../../marcelozerillo/SP80-3280_MMZ.gff3

htseq-count -f sam -r pos -t exon -i Name= -c L2_ScXa_i_10_ntd_count.csv -n 6 L2_ScXa_i_10_ntd.sam ../../marcelozerillo/SP80-3280_MMZ.gff3

htseq-count -f sam -r pos -t gene -i Name= -c L2_ScXa_i_10_ntd_count.csv -n 6 L2_ScXa_i_10_ntd.sam ../../marcelozerillo/SP80-3280_MMZ.gff3
# Funcionou 18/03/2022
htseq-count -f sam -r pos -t gene -i Name -c L2_ScXa_i_10_ntd_count.csv -n 6 L2_ScXa_i_10_ntd.sam ../../marcelozerillo/SP80-3280_MMZ.gff3

