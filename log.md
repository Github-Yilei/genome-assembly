rm commands.list  EVidenceModeler_combined.gff3 nohup.out
nohup sh run.sh &

c001
c002
c003
du ~/project/Geome_assembel/Papeda/05-GenomeAnnotation/Monoploid/Genomethreader/all_pep.blast	147456
du ~/project/Geome_assembel/Papeda/05-GenomeAnnotation/Hap1/03-Genomethreader/all_pep.blast	263936
du ~/project/Geome_assembel/Papeda/05-GenomeAnnotation/Hap2/03-Genomethreader/all_pep.blast	262144



unmasked： EVM合并、从头预测（gene数目和长度都增加了）
augustsus：unmasked.fa的可以使最终结果的 gene数目和长度增加。
EVM合并： masked.fa 可以使最终结果的 gene数目和长度增加。


Gmes: unmasked.fa 结果中，Monoploid gene数目和长度都增加了， 但hap 中的剧烈减少。softmasked 长度和数量表现不如N maksed。

softmasked：蛋白质预测、转录本预测https://github.com/oushujun/EDTA/issues/166

AUGUSTUS： parameter softmasking=1 is slightly more accurate than hard masking (with N's), which looses information. 



EVM： masked  genome获得的基因数目和平均长度均高

简单重复序列可以以很高的统计学显著性比对到蛋白质的低复杂性片段，从而造成假同源。
复杂重复序列包含真正的蛋白质编码基因，会干扰从头预测过程。比如，gene predictor会将位于物种特异性蛋白质编码基因 内含子部分的转座元件误认为是该基因的一部分外显子。
所以识别、标识基因组重复区域非常重要。



Haplotype-resolved genome assembly provides
insights into evolutionary history of the tea plant
Camellia sinensis
