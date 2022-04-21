# QuickTools
目前分为四个板块。


############
一.Quick Plot
可进行GO气泡图、GO柱状图、火山图的一键绘制，输入数据均为xulab RNA-seq 流程产出数据。
绘制步骤：
1.点击打开文件Open按钮，选择输入文件
（气泡图和柱状图需为goea产出的table-de.up-Parent-Child-Union-Benjamini-Hochberg.txt文件，
火山图需为rnadiff产出的rnaseq.sum.tsv）；
2.点击保存路径save按钮选择要保存的文件夹；
3.根据需求选择功能；
4.输入图片的长度和宽度（可选，如不输入，默认高8宽6）；
5.点击start运行。

二、Quick Data Processing
进行一些简单的基因组转录组数据处理
fasta功能：
根据基因id在fasta文件中批量提取序列。打开文件1，选择基因id列表文件，打开文件2选择目标序列fasta文件，
点击start运行。
split功能：
将很大的基因序列fasta文件拆分成若干个小文件
打开文件1，选择大基因序列fasta文件，在辅助参数一栏填写你想要在每个小文件里放入多少条序列
（如，你想将一个包含1000000条fasta序列文件拆分为每个小文件里包含5000条序列，就在辅助参数一栏填写5000）
点击start运行。
TPM、CPM、FPKM计算：
打开文件1（需为xulab rnaseq featurecounts中产出的rnaseq.raw.tsv）
点击start运行

三、Quick Find
根据一列数据，在一个包含这一列数据的大文件中，提取含有该列数据的全部行数据。
也可以取两个单列数据文件的交集部分。
打开单列文件；
打开多列文件；
选择保存路径；
选择分隔符（默认制表符）
填单列对应多列的列索引数（可不填，默认为第一列）
点击start运行

四、Quick Extract
该功能可根据基因ID批量提取基因序列，包括CDS序列，gene序列（含内含子），蛋白序列，
并可以在提取基因序列时，选择上下游序列的长度。
