import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputfile", help="inputfile")
parser.add_argument("-s", "--size", help="size")
parser.add_argument("-o", "--outputfile", help="outputfile")
args = parser.parse_args()

inputDirectory = os.path.abspath(args.inputfile)
outputDirectory = os.path.abspath(args.outputfile)
Size=os.path.abspath(args.size)

i=open(args.inputfile,'r')
alli=i.readlines()
s=open(args.size,'r')
alls=s.readlines()

w=open(args.outputfile,'w')
w.write("chrome"+"\t"+"exon"+"\t"+"intron"+"\t"+"intergenic region"+"\n")

for lines in alls:
    exon=set()
    for linei in alli:
        if linei.split("\t")[2] == lines.split("\t")[0]:
            for i in range(eval(linei.split("\t")[8])):
                exon |= set(range(eval(linei.split("\t")[9].split(",")[i]),eval(linei.split("\t")[10].split(",")[i])))

                # set函数不仅仅用来去重，还用来计算长度（range方法展开成一个个数字代表每个碱基）
                # range eval 函数是用来执行一个字符串表达式#
                #|= reads the same way as +=.

    perexon=float(len(exon))/eval(lines.split("\t")[1])
    #外显子覆盖度是 外显子区域总长/基因组总长
    inter=set()
    for linei in alli:
        if linei.split("\t")[2] == lines.split("\t")[0]:
            inter |= set(range(eval(linei.split("\t")[4]),eval(linei.split("\t")[5])))
    perinter=1-float(len(inter))/eval(lines.split("\t")[1])

    if perexon == 0.0:
        perintron=0.0
        perinter=0.0
    else:
        perintron=1-perexon-perinter

    w.write(lines.split("\t")[0]+"\t"+str(perexon)+"\t"+str(perintron)+"\t"+str(perinter)+"\n")
