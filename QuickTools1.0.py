import sys
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QApplication, QWidget, QTabWidget, QLabel, QLineEdit, QDateEdit,\
                            QComboBox, QTextEdit, QGridLayout, QPushButton,QMessageBox, QVBoxLayout,QFileDialog
import pandas as pd
import numpy as np
from plotnine import *
import re
import json

class Demo(QTabWidget):
    def __init__(self):
        super(Demo, self).__init__()

        self.setWindowTitle('Quick Tools 1.3.2')  #set windows title
        self.setMinimumSize(600, 600)  #set windows size
        self.path = ''                                  # 2
        self.file=''
        self.file2=''
        self.model=''

#fasta function region property
        self.func=''
        self.assist=''




        self.tab1 = QWidget()                   # 1
        self.tab2 = QWidget()
        self.tab3 = QWidget()
        self.tab4 = QWidget()

        self.tab1_init()                        # 2
        self.tab2_init()
        self.tab3_init()
        self.tab4_init()

        self.addTab(self.tab1, 'Quick Plot')    # 3
        self.addTab(self.tab2, 'Quick Data Processing')
        self.addTab(self.tab3, 'Quick Find')
        self.addTab(self.tab4, 'Quick Extract')
    

        self.currentChanged.connect(lambda: print(self.currentIndex()))    # 4

        
    def tab1_init(self):
        open_label = QLabel('打开文件', self.tab1)
        save_label = QLabel('保存路径', self.tab1)
        model_label=QLabel('选择功能', self.tab1)
        start_label = QLabel('启动按钮:', self.tab1)
        height_label=QLabel('图片长度:', self.tab1)
        width_label=QLabel('图片宽度:', self.tab1)

        button_1=QPushButton('Open', self.tab1)
        button_1.clicked.connect(self.open_file_func)
        button_2=QPushButton('Save', self.tab1)
        button_2.clicked.connect(self.save_file_func)
        button_3 = QPushButton('start', self.tab1)
        button_3.clicked.connect(self.get_start1)


        items = ["GO气泡图","GO柱形图","火山图"]
        global model_combo
        model_combo = QComboBox(self.tab1)
        model_combo.addItems(items)
        global height_line
        height_line = QLineEdit(self.tab1)
        global width_line
        width_line = QLineEdit(self.tab1)
        #self.width=width_line.text()

        g_layout = QGridLayout()
        g_layout.addWidget(open_label, 0, 0, 1, 1)
        g_layout.addWidget(button_1, 0, 1, 1, 1)
        g_layout.addWidget(save_label, 1, 0, 1, 1)
        g_layout.addWidget(button_2, 1, 1, 1, 1)
        g_layout.addWidget(model_label, 2, 0, 1, 1)
        g_layout.addWidget(model_combo, 2, 1, 1, 1)
        g_layout.addWidget(height_label, 3, 0, 1, 1)
        g_layout.addWidget(height_line, 3, 1, 1, 1)
        g_layout.addWidget(width_label, 4, 0, 1, 1)
        g_layout.addWidget(width_line, 4, 1, 1, 1)
        g_layout.addWidget(start_label, 5, 0, 1, 1)
        g_layout.addWidget(button_3, 5, 1, 1, 1)

        self.tab1.setLayout(g_layout)

    def tab2_init(self):
        open_label = QLabel('打开文件1', self.tab2)
        open_label2 = QLabel('打开文件2', self.tab2)
        save_label = QLabel('保存路径', self.tab2)
        func_label=QLabel('选择功能', self.tab2)
        assist_label=QLabel('辅助参数', self.tab2)
        start_label = QLabel('启动按钮:', self.tab2)

        button_open1=QPushButton('Open', self.tab2)
        button_open1.clicked.connect(self.open_file_func)
        button_open2=QPushButton('Open', self.tab2)
        button_open2.clicked.connect(self.open_file2_func)
        button_save=QPushButton('Save', self.tab2)
        button_save.clicked.connect(self.save_file_func)
        button_start = QPushButton('start', self.tab2)
        button_start.clicked.connect(self.get_start2)

        items = ["fasta","split","CPM","TPM","FPKM"]
        global func_combo
        func_combo = QComboBox(self.tab2)
        func_combo.addItems(items)
        global Assist_line
        Assist_line=QLineEdit(self.tab2)

        g_layout = QGridLayout()
        g_layout.addWidget(open_label, 0, 0, 1, 1)
        g_layout.addWidget(button_open1, 0, 1, 1, 1)
        g_layout.addWidget(open_label2, 1, 0, 1, 1)
        g_layout.addWidget(button_open2, 1, 1, 1, 1)
        g_layout.addWidget(save_label, 2, 0, 1, 1)
        g_layout.addWidget(button_save, 2, 1, 1, 1)
        g_layout.addWidget(func_label, 3, 0, 1, 1)
        g_layout.addWidget(func_combo, 3, 1, 1, 1)
        g_layout.addWidget(assist_label, 4, 0, 1, 1)
        g_layout.addWidget(Assist_line, 4, 1, 1, 1)
        g_layout.addWidget(start_label, 5, 0, 1, 1)
        g_layout.addWidget(button_start, 5, 1, 1, 1)

        self.tab2.setLayout(g_layout)

#####tab3
    def tab3_init(self):
        open_label = QLabel('打开单列文件', self.tab3)
        open_label2 = QLabel('打开多列文件', self.tab3)
        save_label = QLabel('保存路径', self.tab3)
        sep_label=QLabel('选择分隔符', self.tab3)
        col_label=QLabel('选择列数', self.tab3)
        start_label = QLabel('启动按钮:', self.tab3)

        button_open1=QPushButton('Open', self.tab3)
        button_open1.clicked.connect(self.open_file_func)
        button_open2=QPushButton('Open', self.tab3)
        button_open2.clicked.connect(self.open_file2_func)
        button_save=QPushButton('Save', self.tab3)
        button_save.clicked.connect(self.save_file_func)
        button_start = QPushButton('start', self.tab3)
        button_start.clicked.connect(self.get_start3)

        items = ["制表符","逗号","空格"]
        global sep_line
        sep_line=QComboBox(self.tab3)
        sep_line.addItems(items)
        global col_line
        col_line=QLineEdit(self.tab3)

        g_layout = QGridLayout()
        g_layout.addWidget(open_label, 0, 0, 1, 1)
        g_layout.addWidget(button_open1, 0, 1, 1, 1)
        g_layout.addWidget(open_label2, 1, 0, 1, 1)
        g_layout.addWidget(button_open2, 1, 1, 1, 1)
        g_layout.addWidget(save_label, 2, 0, 1, 1)
        g_layout.addWidget(button_save, 2, 1, 1, 1)
        g_layout.addWidget(sep_label, 3, 0, 1, 1)
        g_layout.addWidget(sep_line, 3, 1, 1, 1)
        g_layout.addWidget(col_label, 4, 0, 1, 1)
        g_layout.addWidget(col_line, 4, 1, 1, 1)
        g_layout.addWidget(start_label, 5, 0, 1, 1)
        g_layout.addWidget(button_start, 5, 1, 1, 1)

        self.tab3.setLayout(g_layout)

######tab4
    def tab4_init(self):
        open_label = QLabel('打开文件', self.tab4)
        save_label = QLabel('保存路径', self.tab4)
        genome_label=QLabel('选择基因组',self.tab4)
        ex_part_label=QLabel('选择模式', self.tab4)
        up_label=QLabel('选择上游序列', self.tab4)
        down_label=QLabel('选择下游序列', self.tab4)
        start_label = QLabel('启动按钮:', self.tab4)

        button_open1=QPushButton('Open', self.tab4)
        button_open1.clicked.connect(self.open_file_func)
        button_save=QPushButton('Save', self.tab4)
        button_save.clicked.connect(self.save_file_func)
        button_start = QPushButton('start', self.tab4)
        button_start.clicked.connect(self.get_start_ex)

        genos=["FGRRES",'Neurospora_crassa']
        global genome_combo
        genome_combo=QComboBox(self.tab4)
        genome_combo.addItems(genos)
        items = ["CDS","gene","protein"]
        global ex_part_combo
        ex_part_combo=QComboBox(self.tab4)
        ex_part_combo.addItems(items)
        global up_line
        up_line=QLineEdit(self.tab4)
        global down_line
        down_line=QLineEdit(self.tab4)


        g_layout = QGridLayout()
        g_layout.addWidget(open_label, 0, 0, 1, 1)
        g_layout.addWidget(button_open1, 0, 1, 1, 1)
        g_layout.addWidget(save_label, 1, 0, 1, 1)
        g_layout.addWidget(button_save, 1, 1, 1, 1)
        g_layout.addWidget(genome_label, 2, 0, 1, 1)
        g_layout.addWidget(genome_combo, 2, 1, 1, 1)
        g_layout.addWidget(ex_part_label, 3, 0, 1, 1)
        g_layout.addWidget(ex_part_combo, 3, 1, 1, 1)
        g_layout.addWidget(up_label, 4, 0, 1, 1)
        g_layout.addWidget(up_line, 4, 1, 1, 1)
        g_layout.addWidget(down_label, 5, 0, 1, 1)
        g_layout.addWidget(down_line, 5, 1, 1, 1)
        g_layout.addWidget(start_label, 6, 0, 1, 1)
        g_layout.addWidget(button_start, 6, 1, 1, 1)

        self.tab4.setLayout(g_layout)





#######plot function region 
    def open_file_func(self):                           # 5
        file, _ = QFileDialog.getOpenFileName(self, 'Open File', './', 'Files (*)')
        if file:
            self.file=file

    def open_file2_func(self):                           # 5
        file, _ = QFileDialog.getOpenFileName(self, 'Open File', './', 'Files (*)')
        if file:
            self.file2=file

    def save_file_func(self):
        file=QFileDialog.getExistingDirectory(self, 'Select FileDirectory', './')
        if file:
            self.path=file



    def get_model(self):
        self.model=model_combo.currentText()  #combox.currentText()读取当前选项的值变为字符串
        return self.model

    def get_height(self):
        height=height_line.text()
        return height
    
    def get_width(self):
        width=width_line.text()  #.text 返回 QLineEdit当前值 但得单独写函数
        return width

    def plot_bubble(self):
        file=self.file
        df=pd.read_table(file,sep='\t')
        df=df[df['p.adjusted']<0.05]
        df=df.sort_values(by=['name.space','p'])
        df['name.space']=df['name.space'].astype('category')
        gofold=(df['Study.term']/df['Study.total'])/(df['Pop.term']/df['Pop.total'])
        gof=df['p.adjusted'].map(lambda x:-np.log10(x))
        p1=ggplot(df,aes(x=gofold,y='name',color='Study.term',size=gof))+geom_point()+labs(color='Counts',alpha="Number of enriched genes",size="−log10(FDR)",x="Enrichment factor (fold)",y="",title="")+scale_size_continuous(range=(4,8))+theme_bw()+theme_classic()+ theme(figure_size=(8, 6))+ ylim(df.iloc[:,12])

        raw_height=self.get_height()
        if raw_height != '':
            height=raw_height
        else:
            height=8
        raw_width=self.get_width()
        if raw_width != '':
            width=raw_width
        else:
            width=6
        p1.save(r'%s\bubble.pdf'%(self.path), height=int(height), width=int(width))


    def GO_barplot(self):
        file=self.file
        df=pd.read_table(file,sep='\t')
        df['p']=df['p'].map(lambda x:-np.log10(x))
        df['p']=df['p'].round(2)
        df=df[df['p.adjusted']<0.05]
        #df=df.sort_values(by='p')
        df['name.space'].replace({'biological_process':'BP','molecular_function':'MF','cellular_component':'CC'},inplace=True)
        df=df.sort_values(by=['name.space','p'])
        df['p']=df['p'].astype('category')
        df['name.space']=df['name.space'].astype('category')
        CPCOLS =('#8DA1CB','#FD8D62','#66C3A5')
        p1 = ggplot(df,aes(x='name',y='p',fill='name.space'))+geom_bar(stat = 'identity',width = 0.3)+coord_flip() +scale_fill_manual(values = CPCOLS)+theme_classic()+ ylab('-log(P value)')+labs(fill='Category')+ xlim(df.iloc[:,12])
        raw_height=self.get_height()
        if raw_height != '':
            height=raw_height
        else:
            height=8
        raw_width=self.get_width()
        if raw_width != '':
            width=raw_width
        else:
            width=6
        p1.save(f'{self.path}\\barplot.pdf', height=int(height), width=int(width))

        
    def volcano_plot(self):
        file=self.file
        df=pd.read_table(file,sep='\t',index_col=0)
        index_list=[i for i in range(len(df.index))]
        f=lambda x:'up' if df.iloc[x,3]<0.05 and df.iloc[x,0]>=1 else('Down' if df.iloc[x,3]<0.05 and df.iloc[x,0]<=-1 else('Stable'))
        df['change']=list(map(f,index_list)) 
        df['FDR']=df['FDR'].map(lambda x:-np.log10(x))
        p1=ggplot(df,aes(x='logFC',y='FDR',colour='change'))+geom_point(alpha=1, size=2)+scale_color_manual(values=("red", "#00B2FF","orange"))+geom_vline(xintercept=[-1,1],color="#990000",size=0.8,linetype="dashed")+geom_hline(yintercept = -np.log10(0.05),color="#990000",size=0.8,linetype="dashed")+labs(x="log2(fold change)",y="-log10 (FDR)")+theme_bw()
        raw_height=self.get_height()
        if raw_height != '':
            height=raw_height
        else:
            height=8
        raw_width=self.get_width()
        if raw_width != '':
            width=raw_width
        else:
            width=6
        p1.save(f'{self.path}\\volcanoPlot.pdf', height=int(height), width=int(width))





######fasta region
    def get_func(self):
        self.func=func_combo.currentText()  #combox.currentText()读取当前选项的值变为字符串
        return self.func

    def get_assist(self):
        self.assist=Assist_line.text()  #.text 返回 QLineEdit当前值 但得单独写函数
        return self.assist


    def fasta(self):
        file1=self.file
        file2=self.file2
        fr=open(file2,'r')
        fl=open(file1,'r')
        fw=open(f'{self.path}\\output.fasta','w')
        dict={}
        for line in fr:
            if line.startswith('>'):
                id_is = line.replace('>','').split()[0]
                dict[id_is]=''
            else:
                dict[id_is]+=line.replace('\n','')
        fr.close()
        for id in fl:
            id=id.strip()
            for key in dict.keys():
                m=re.search(id,key)
                if m is not None:
                    fw.write('>'+id)
                    fw.write('\n')
                    fw.write(dict[key])
                    fw.write('\n')
        fl.close()
        fw.close()

    def split(self):
        input_file=self.file
        prefix='output'
        num_m=self.get_assist()
        X=int(num_m)
        FA_in_file = open(input_file, "r")
        fa_Info = []
        fa_Seq = []
        fa_Num = -1
        for Y in FA_in_file.readlines():
            Y = Y.rstrip()
            if Y[0] == ">":
                fa_Info.append(Y)
                fa_Num = fa_Num + 1
                fa_Seq.append("")
            else:
                fa_Seq[fa_Num] = fa_Seq[fa_Num] + Y
        file_Num = (fa_Num + 1)//X + 1
        for i in range(file_Num):
            exec(prefix + str(i + 1) + ' = open("' +self.path+'\\'+prefix + str(i + 1) + '.fasta"' + ', "w")')
            start = i * X
            end = (i + 1) * X
            if end > fa_Num + 1:
                end = fa_Num + 1
            for j in range(start, end, 1):
                exec(prefix + str(i + 1) + '.write(fa_Info[j] + "\\n")')
                while len(fa_Seq[j]) > 60:
                    exec(prefix + str(i + 1) + '.write(fa_Seq[j][:60] + "\\n")')
                    fa_Seq[j] = fa_Seq[j][60:]
                else:
                    exec(prefix + str(i + 1) + '.write(fa_Seq[j] + "\\n")')
            exec(prefix + str(i + 1) + '.close()')
        FA_in_file.close()





    def cpm(self):
        file=self.file
        foldchange=pd.read_table(file,skiprows=1)
        rownum=foldchange.shape[1]
        countdata=foldchange.iloc[:,6:rownum]
        colsumn=countdata.apply(lambda x:x.sum(),axis=0)#axis = 0 按行计算,得到列的性质。axis = 1 按列计算,得到行的性质。
        cpm=countdata/colsumn*1000000
        col=foldchange.iloc[:,0]
        cpm=pd.concat([col,cpm],axis=1)
        cpm.to_csv(f'{self.path}\\cpm.csv',sep=',',index=False)

    def tpm(self):
        file=self.file
        foldchange=pd.read_table(file,skiprows=1)
        rownum=foldchange.shape[1]
        countdata=foldchange.iloc[:,6:rownum]
        kb=foldchange.iloc[:,5]/1000
        rpk=countdata.div(kb,axis=0)
        tpm=rpk/rpk.sum()*1000000
        col=foldchange.iloc[:,0]
        tpm=pd.concat([col,tpm],axis=1)
        tpm.to_csv(f'{self.path}\\tpm.csv',sep=',',index=False)

    def fpkm(self):
        file=self.file
        foldchange=pd.read_table(file,skiprows=1)
        rownum=foldchange.shape[1]
        countdata=foldchange.iloc[:,6:rownum]
        kb=foldchange.iloc[:,5]/1000
        rpk=countdata.div(kb,axis=0)
        fpkm=rpk/countdata.sum()*1000000
        col=foldchange.iloc[:,0]
        fpkm=pd.concat([col,fpkm],axis=1)
        fpkm.to_csv(f'{self.path}\\fpkm.csv',sep=',',index=False)
#####tab3 func
    def get_sep_col(self):
        sep=sep_line.currentText()
        col=col_line.text()
        return sep, col


    def find(self):
        sep,col=self.get_sep_col()
        if sep=="制表符":
            sep_='\t'
        elif sep=="逗号":
            sep_=","
        elif sep=="空格":
            sep_=' '
        F=0
        if col !='':
            F=int(col)-1
        print(F)
        file1=self.file
        file2=self.file2
        idlist=open(file1,'r')
        findfile=open(file2,'r')
        format_=file2.strip().split('.')[-1]
        outfile=open(f'{self.path}\\target.{format_}','w')
        mydic={}
        for m in findfile:
            mchange=m.strip().split(sep_)[F]
            mydic[mchange]=m
        for i in idlist:
            i=i.strip()
            for key in mydic.keys():
                if key==i:
                    outfile.write(mydic[key])
        idlist.close()
        findfile.close()
        outfile.close()
#############################
#############################
    def get_genome(self):
        genome=genome_combo.currentText()
        return genome


    def readgff(self):
        genome=self.get_genome()
        if genome=='FGRRES':
            with open('.\\jsonFile\\fgrres_gff.json','r') as gff:
                cds=json.load(gff)
        if genome=='Neurospora_crassa':
            with open('.\\jsonFile\\ncgtf.json','r') as gff:
                cds=json.load(gff)
        return cds
    def readfa(self):
        genome=self.get_genome()
        if genome=='FGRRES':
            with open('.\\jsonFile\\fgrr_fa.json','r') as fa:
                chr=json.load(fa)
        if genome=='Neurospora_crassa':
            with open('.\\jsonFile\\ncfa.json','r') as gff:
                chr=json.load(gff)
        return chr

#找到cds的开头和结尾位置
    def cdsM(self):
        cds=self.readgff()
        diccdsM={}
        for key in cds.keys():
            diccdsM[key]=[cds[key][0][0],cds[key][0][1],cds[key][-1][2],cds[key][-1][3]]
        return diccdsM

    def exCDS(self):
        
        
        position=self.readgff()
        cdsdic={}
        with open(self.file,'r',encoding='UTF-8') as genelist:
            for line in genelist.readlines():
                l=line.strip()
                cdsdic[l]=''
                
                if position[l][0][-1]=='+':
                    for i in position[l]:
                        cdsdic[l]+=str(self.seq(i))
                if position[l][0][-1]=='-':
                    for i in position[l][::-1]:
                        cdsdic[l]+=str(self.seq(i))

        return cdsdic

    def output_CDS(self):
        cdsdic=self.exCDS()
        with open(f'{self.path}\\cds.fasta','w') as outfile:
            out=cdsdic
            for key in out.keys():
                outfile.write('>'+key+'\n'+out[key]+'\n')
                
    def exGENE(self):
        genedic={}
        cdsM=self.cdsM()
        with open(self.file,'r',encoding='UTF-8') as genelist:
            for line in genelist.readlines():
                l=line.strip()
                genedic[l]=str(self.seq(cdsM[l]))
        with open(f'{self.path}\\gene.fasta','w') as outfile:
            out=genedic
            for key in out.keys():
                outfile.write('>'+key+'\n'+out[key]+'\n')
        return
    

    def seq(self,position):
        raw_up,raw_down=self.up_down()
        if raw_up!='':
            up=raw_up
        else:
            up=0
        if raw_down!='':
            down=raw_down
        else:
            down=0
        al=self.readfa()
        hash={'A':'T','T':'A','C':'G','G':'C','N':'N'}
        seqid,start,end,strand=position
        if strand=='+':
            dna=al[seqid][int(start)-1-int(up):int(end)+int(down)].upper()
        
        if strand=='-':
            dnar=al[seqid][int(start)-1-int(down):int(end)+int(up)][::-1].upper()
            dna=''.join([hash[i] for i in dnar])
        return dna
#根据cds提取蛋白序列
    def exPRO(self):
        cds=self.exCDS()
        aa_dict = {
    'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I','ATA':'I', 
    'ATG':'M', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'CCT':'P', 'CCC':'P',
    'CCA':'P', 'CCG':'P', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A','TAT':'Y',
    'TAC':'Y',  'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 
    'GAA':'E', 'GAG':'E', 'TGT':'C', 'TGC':'C', 'TGG':'W', 'CGT':'R', 'CGC':'R','CGA':'R','CGG':'R', 'AGT':'S', 'AGC':'S', 
    'AGA':'R', 'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'TGA':'*', 'TAA':'*', 'TAG':'*',
    }
        aa_pro={}    
        for k in cds.keys():
            aa_seq=[]
            for i in range(0,len(cds[k]),3):
                codens=cds[k][i:i+3]
                if aa_dict.get(codens):#确定codens是否是氨基酸字典里的键,if后，非0非空为真，0或None为假
                    aa_seq.append(aa_dict[codens])
            aa=''.join(aa_seq).replace('*','')
            aa_pro[k]=aa
        with open(f'{self.path}\\pep.fasta','w') as outfile:
            out=aa_pro
            for key in out.keys():
                outfile.write('>'+key+'\n'+out[key]+'\n')
        return

    def get_ex_part(self):
        ex=ex_part_combo.currentText()
        return ex
    
    def up_down(self):
        up=up_line.text()
        down=down_line.text()
        return up,down
        
    def get_start_ex(self):
        ex_part=self.get_ex_part()

        if ex_part=="CDS":
            self.output_CDS()
        elif ex_part=="gene":
            self.exGENE()
        elif ex_part=="protein":
            self.exPRO()
















#########################

    def get_start1(self):
        model=self.get_model()

        if model=="GO气泡图":
            self.plot_bubble()
        elif model=="GO柱形图":
            self.GO_barplot()
        elif model=="火山图":
            self.volcano_plot()   
    def get_start2(self):
        func=self.get_func()
        if func=="fasta":
            self.fasta()
        elif func=="split":
            self.split()
        elif func=="CPM":
            self.cpm()
        elif func=="TPM":
            self.tpm()
        elif func=="FPKM":
            self.fpkm()

    def get_start3(self):
        self.find()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = Demo()
    demo.show()
    sys.exit(app.exec_())
