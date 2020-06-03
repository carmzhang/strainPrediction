import xlwt
import os

from optparse import OptionParser
import subprocess
import sys
import difflib
import xlsxwriter
import re
import math
path="D:\\Duke\\Duke_201904\\Resfinder_0.8_0.6_20190426"
path3="D:\\Duke\\Duke_201904"

def getwork():
    files=os.listdir(path)
    wb=xlsxwriter.Workbook((path3+"\\"+'Resfinder_results20190426_640.xlsx'))
    ws1 = wb.add_worksheet('Resfinder')
    dic_tmp={}
    dic_ID={}
    results=[]
    idlist={}
    nb=1
    ws1.write (0,0,'ID')
    listres=["Sulphonamide","MLS - Macrolide, Lincosamide and Streptogramin B","Fusidic Acid","Nitroimidazole","Aminoglycoside","Oxazolidinone","Beta-lactam","Fosfomycin","Rifampicin","Colistin","Tetracycline","Glycopeptide","Phenicol","Trimethoprim","Fluoroquinolone"]
    for g in listres:
        dic_ID[g]=nb
        ws1.write(0,nb,g)
        nb = nb+1
    for i in range(len(files)):
        workfile = path+"\\"+files[i]
        file=files[i]
        file_name=file.split('.')
        file_ID=file_name[0].split('_')
        i += 1
        ws1.write(i,0,str(file_ID[0]))
        idlist[file_ID[0]]=i
        infile = open(workfile, "r").readlines()
        number=0
        for gl in listres:
            dic_tmp[gl]=[]
            for m in range(0,len(infile)):
                line1=infile[m].strip()
                if line1.startswith('Resistance'):
                    naline=infile[m-1].strip()
                    if naline==gl:
                        while not infile[m+1].startswith('\n'):
                            idco=infile[m+1].split("\t")
                            m += 1
                            identity=idco[1]
                            coverage=idco[3]
                            genelength=idco[2].split("/")[1]
                            results=[str(idco[0]),str(identity),str(coverage),str(genelength),'|']
                            dic_tmp[gl].extend(results)
                            #print(file_ID[0],gl,dic_tmp[gl],"!!!",m,idco)
                            results=[]
                            
        for ider in dic_tmp.keys():
            tmp=(' ').join(dic_tmp[ider])
            ws1.write(i,dic_ID[ider],tmp)
    wb.close()
    
if __name__=='__main__':
    getwork()
